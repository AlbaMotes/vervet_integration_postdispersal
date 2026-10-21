library(dplyr)
library(purrr)
library(tibble)
library(brms)
library(posterior)


group_presence_matrices<-readRDS("/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/Pooja_chapter_matrank/Data/group_presence_matrices.rds")
emigrants_full <- read.csv(
  "/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/Pooja_chapter_matrank/Data/emigrants_full_merged.csv"
)

life_hist <- read.csv(
  "/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/IVP_data/life_history.csv"
)
# ── Helper: Codes present (==1) in a given group on a given date ───────
# Falls back to the nearest available date if the exact date isn't a row
# in the presence matrix (matrices are built daily from group start, but
# defensive fallback avoids silent NA on any edge-case gaps).
get_present_codes <- function(group_name, target_date, group_presence_matrices) {
  if (is.na(group_name) || is.na(target_date) ||
      !group_name %in% names(group_presence_matrices)) {
    return(character(0))
  }
  pm  <- group_presence_matrices[[group_name]]
  row <- pm %>% filter(Date == as.Date(target_date))
  if (nrow(row) == 0) {
    row <- pm %>%
      mutate(diff = abs(as.numeric(Date - as.Date(target_date)))) %>%
      arrange(diff) %>% slice(1) %>% select(-diff)
  }
  if (nrow(row) == 0) return(character(0))
  names(row)[-1][as.logical(unlist(row[1, -1]) == 1)]
}

check_convergence <- function(fit, rhat_threshold = 1.01, ess_threshold = 400) {
  summ <- summarise_draws(as_draws_df(fit))
  bad_rhat <- summ %>% filter(rhat > rhat_threshold)
  bad_ess  <- summ %>% filter(ess_bulk < ess_threshold | ess_tail < ess_threshold)
  ok <- nrow(bad_rhat) == 0 && nrow(bad_ess) == 0
  if (!ok) {
    cat("  CONVERGENCE WARNING:\n")
    if (nrow(bad_rhat) > 0) print(bad_rhat %>% select(variable, rhat))
    if (nrow(bad_ess)  > 0) print(bad_ess  %>% select(variable, ess_bulk, ess_tail))
  }
  list(ok = ok, max_rhat = max(summ$rhat, na.rm = TRUE),
       min_ess = min(c(summ$ess_bulk, summ$ess_tail), na.rm = TRUE))
}

library(tibble)

summarize_posterior <- function(estimates) {
  tibble(
    estimate_mean   = mean(estimates),
    estimate_median = median(estimates),
    network_q025    = quantile(estimates, 0.025, names = FALSE),
    network_q975    = quantile(estimates, 0.975, names = FALSE),
    prop_positive   = mean(estimates > 0),
    prop_negative   = mean(estimates < 0)
  )
}
# ── Helper: Codes present at any point across a date range in a group ──
get_present_codes_range <- function(group_name, start_date, end_date, group_presence_matrices) {
  if (is.na(group_name) || is.na(start_date) || is.na(end_date) ||
      !group_name %in% names(group_presence_matrices)) {
    return(character(0))
  }
  pm <- group_presence_matrices[[group_name]] %>%
    filter(Date >= as.Date(start_date), Date <= as.Date(end_date))
  if (nrow(pm) == 0) return(character(0))
  present_cols <- names(pm)[-1][colSums(pm[, -1, drop = FALSE] == 1, na.rm = TRUE) > 0]
  present_cols
}

sex_lookup <- life_hist %>%
  mutate(Code = tolower(trimws(Code))) %>%
  select(Code, Sex) %>%
  distinct()

# ══════════════════════════════════════════════════════════════════════
# Compute all three variables in one pass over emigrants_full
# ══════════════════════════════════════════════════════════════════════
demo_vars <- pmap_dfr(
  list(
    Code             = emigrants_full$Code,
    BirthGp          = emigrants_full$BirthGp,
    ImmigrationGp1   = emigrants_full$ImmigrationGp1,
    DOB              = as.Date(emigrants_full$DOB),
    DateImmigration1 = as.Date(emigrants_full$DateImmigration1),
    DepartureNatalGp = as.Date(emigrants_full$DepartureNatalGp),
    juvenile_period_end = as.Date(emigrants_full$juvenile_period_end)
  ),
  function(Code, BirthGp, ImmigrationGp1, DOB, DateImmigration1,
           DepartureNatalGp, juvenile_period_end) {
    
    dispersal_date <- if (!is.na(DateImmigration1)) DateImmigration1 else DepartureNatalGp
    
    # ── 1. Age at dispersal, in months ─────────────────────────────────
    age_months <- if (!is.na(dispersal_date) && !is.na(DOB)) {
      as.numeric(dispersal_date - DOB) / (365.25 / 12)
    } else {
      NA_real_
    }
    
    # ── 2. Sex ratio of destination group on dispersal date ────────────
    dest_present <- get_present_codes(ImmigrationGp1, dispersal_date, group_presence_matrices)
    dest_sexes   <- sex_lookup$Sex[match(dest_present, sex_lookup$Code)]
    n_males_dest   <- sum(dest_sexes == "M", na.rm = TRUE)
    n_females_dest <- sum(dest_sexes == "F", na.rm = TRUE)
    sex_ratio_dest <- if (length(dest_present) == 0) {
      NA_real_
    } else if (n_females_dest == 0) {
      NA_real_   # undefined ratio; flagged rather than Inf
    } else {
      n_males_dest / n_females_dest
    }
    
    # ── 3. Known males from natal juvenile period, present in destination ─
    natal_present <- get_present_codes_range(BirthGp, DOB, juvenile_period_end, group_presence_matrices)
    natal_sexes   <- sex_lookup$Sex[match(natal_present, sex_lookup$Code)]
    known_males   <- setdiff(natal_present[natal_sexes == "M"], Code)  # exclude self
    
    n_known_males_in_dest <- if (length(dest_present) == 0) {
      NA_integer_
    } else {
      length(intersect(known_males, dest_present))
    }
    
    tibble(
      Code                   = Code,
      age_at_dispersal_months = age_months,
      dest_n_males           = if (length(dest_present) == 0) NA_integer_ else n_males_dest,
      dest_n_females         = if (length(dest_present) == 0) NA_integer_ else n_females_dest,
      dest_sex_ratio         = sex_ratio_dest,
      n_known_males_in_dest  = n_known_males_in_dest
    )
  }
)

emigrants_full <- emigrants_full %>%
  left_join(demo_vars, by = "Code")

# ── Diagnostics ──────────────────────────────────────────────────────
cat("=== Missingness summary ===\n")
emigrants_full %>%
  summarise(
    n_total                    = n(),
    n_age_at_dispersal         = sum(!is.na(age_at_dispersal_months)),
    n_dest_sex_ratio           = sum(!is.na(dest_sex_ratio)),
    n_dest_sex_ratio_undefined = sum(!is.na(dest_n_males) & dest_n_females == 0),
    n_known_males_in_dest      = sum(!is.na(n_known_males_in_dest))
  ) %>% print()

cat("\n=== age_at_dispersal_months ===\n"); print(summary(emigrants_full$age_at_dispersal_months))
cat("\n=== dest_sex_ratio ===\n");          print(summary(emigrants_full$dest_sex_ratio))
cat("\n=== dest_n_males ===\n");            print(summary(emigrants_full$dest_n_males))
cat("\n=== n_known_males_in_dest ===\n");   print(summary(emigrants_full$n_known_males_in_dest))
cat("\nDistribution of n_known_males_in_dest:\n")
print(table(emigrants_full$n_known_males_in_dest, useNA = "ifany"))

write.csv(emigrants_full,
          file.path(data_dir, "emigrants_full_with_demographic_predictors.csv"),
          row.names = FALSE)


run_demographic_mediation <- function(emigrants_full, mediator_col, mediator_label,
                                      chains = 2, iter = 4000, warmup = 1000,
                                      adapt_delta = 0.95) {
  
  cat("Running:", mediator_label, "\n")
  
  dat <- emigrants_full %>%
    filter(!is.na(mother_avg_elo),
           !is.na(elo_one_year_male_male),
           !is.na(.data[[mediator_col]]),
           !is.na(obs_weight))
  
  if (nrow(dat) < 15) { cat("  Too few males:", nrow(dat), "- skipping\n"); return(NULL) }
  cat("  n =", nrow(dat), "\n")
  
  dat <- dat %>%
    mutate(
      maternal_elo_z = scale(mother_avg_elo)[, 1],
      mediator_z     = scale(.data[[mediator_col]])[, 1],
      male_elo_z     = scale(elo_one_year_male_male)[, 1]
    )
  
  bf_mediator <- bf(mediator_z ~ maternal_elo_z, sigma ~ obs_weight)
  bf_outcome  <- bf(male_elo_z ~ mediator_z + maternal_elo_z, sigma ~ obs_weight)
  
  fit <- brm(
    bf_mediator + bf_outcome + set_rescor(FALSE),
    data = dat, family = gaussian(),
    chains = chains, iter = iter, warmup = warmup,
    seed = 1, refresh = 0, silent = 2,
    control = list(adapt_delta = adapt_delta)
  )
  
  conv <- check_convergence(fit)
  if (!conv$ok) cat("  WARNING: did not converge cleanly (max Rhat =",
                    round(conv$max_rhat, 4), ", min ESS =", round(conv$min_ess, 1), ")\n")
  
  post <- as_draws_df(fit)
  a_post   <- post$b_mediatorz_maternal_elo_z
  b_post   <- post$b_maleeloz_mediator_z
  c_post   <- post$b_maleeloz_maternal_elo_z
  ind_post <- a_post * b_post
  tot_post <- c_post + ind_post   # total effect, same identity as before
  
  results <- bind_rows(
    summarize_posterior(a_post)   %>% mutate(path = "a (maternal rank -> mediator)"),
    summarize_posterior(b_post)   %>% mutate(path = "b (mediator -> male elo)"),
    summarize_posterior(c_post)   %>% mutate(path = "c' (maternal rank -> male elo, direct)"),
    summarize_posterior(ind_post) %>% mutate(path = "indirect (a * b)"),
    summarize_posterior(tot_post) %>% mutate(path = "c (total effect)")
  ) %>%
    mutate(mediator = mediator_label, n = nrow(dat), .before = 1)
  
  cat("\n  === Results for", mediator_label, "(n =", nrow(dat), ") ===\n")
  for (j in 1:nrow(results)) {
    cat(sprintf("  %-45s mean = %6.3f  95%% CI [%6.3f, %6.3f]  pd = %4.1f%%\n",
                results$path[j], results$estimate_mean[j],
                results$network_q025[j], results$network_q975[j],
                max(results$prop_positive[j], results$prop_negative[j]) * 100))
  }
  
  list(summary = results, fit = fit)
}

# ── Run for all three candidate mediators ────────────────────────────
mediator_list <- list(
  age_at_dispersal = "age_at_dispersal_months",
  sex_ratio        = "dest_sex_ratio",
  n_males          = "n_known_males_in_dest"
)

emigrants_full <- emigrants_full %>%
  mutate(
    DateImmigration1 = as.Date(DateImmigration1),
    DepartureNatalGp = as.Date(DepartureNatalGp),
    LastSeen1        = as.Date(LastSeen1),
    start_date_effective = if_else(!is.na(DateImmigration1), DateImmigration1, DepartureNatalGp),
    window_end        = start_date_effective + 365,
    effective_end     = if_else(!is.na(LastSeen1) & LastSeen1 < window_end, LastSeen1, window_end),
    days_present_first_year = pmin(as.numeric(effective_end - start_date_effective), 365),
    obs_weight        = days_present_first_year / 365
  )

all_demo_results <- imap(mediator_list, function(col, label) {
  run_demographic_mediation(emigrants_full, col, label)
})

demo_summary <- bind_rows(map(all_demo_results, "summary"))
print(as.data.frame(demo_summary), row.names = FALSE)

write.csv(demo_summary,
          file.path(data_dir, "demographic_mediation_results_bayesian.csv"),
          row.names = FALSE)