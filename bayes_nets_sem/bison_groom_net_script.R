

library(dplyr)
library(tidyr)
library(igraph)

#The "event" for dyad A→B is the number of days A groomed B, and the "duration" is the total observation days for that pair.
#emigrants<-read.csv("/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/Pooja_chapter_matrank/Data//All_potential_immigrants_Josie_version.csv")
#read grooming_dataset
#grooming_data<-read.csv("/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Pooja_chapter_matrank/Data/grooming_data2010_2025.csv")
#to count number of adlib days per group
#affiliative <- read.csv("/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/IVP_data/adlib/affiliative.csv")
#group presence matrices
#group_presence_matrices <- readRDS("/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/Pooja_chapter_matrank/Data/group_presence_matrices.rds")

group_abbreviations <- c("Baie Dankie" = "BD", "Noha" = "NH", "Ankhase" = "AK", 
                         "Kubu" = "KB", "Lemon Tree" = "LT", "Crossing" = "CR", 
                         "IFamily" = "IF")
affiliative$group_Abbreviation <- group_abbreviations[affiliative$group]

# ══════════════════════════════════════════════════════════════════════
# Helper: compute directed network metrics from Beta posterior samples
# ══════════════════════════════════════════════════════════════════════
compute_directed_metrics_from_posterior <- function(dyadic_data, all_nodes,
                                                    num_draws = 1000,
                                                    prior_a = 0.1, prior_b = 0.1,
                                                    edge_thresh = 0.0001,
                                                    group_size = NULL) {
  
  n_dyads <- nrow(dyadic_data)
  n_nodes <- length(all_nodes)
  
  draw_mat <- matrix(NA, nrow = n_dyads, ncol = num_draws)
  for (d in 1:n_dyads) {
    draw_mat[d, ] <- rbeta(num_draws,
                           prior_a + dyadic_data$event[d],
                           prior_b + dyadic_data$duration[d] - dyadic_data$event[d])
  }
  
  in_strength_mat  <- matrix(NA, nrow = num_draws, ncol = n_nodes, dimnames = list(NULL, all_nodes))
  out_strength_mat <- matrix(NA, nrow = num_draws, ncol = n_nodes, dimnames = list(NULL, all_nodes))
  in_degree_mat    <- matrix(NA, nrow = num_draws, ncol = n_nodes, dimnames = list(NULL, all_nodes))
  out_degree_mat   <- matrix(NA, nrow = num_draws, ncol = n_nodes, dimnames = list(NULL, all_nodes))
  eigen_mat        <- matrix(NA, nrow = num_draws, ncol = n_nodes, dimnames = list(NULL, all_nodes))
  
  from_node <- as.character(dyadic_data$groomer)
  to_node   <- as.character(dyadic_data$groomee)
  
  for (k in 1:num_draws) {
    
    adj <- matrix(0, n_nodes, n_nodes, dimnames = list(all_nodes, all_nodes))
    for (r in 1:n_dyads) {
      w <- draw_mat[r, k]
      if (w < edge_thresh) w <- 0
      adj[from_node[r], to_node[r]] <- w
    }
    
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "directed", weighted = TRUE)
    
    in_strength_mat[k, ]  <- igraph::strength(g, mode = "in")
    out_strength_mat[k, ] <- igraph::strength(g, mode = "out")
    in_degree_mat[k, ]    <- igraph::degree(g, mode = "in")  / (group_size - 1)
    out_degree_mat[k, ]   <- igraph::degree(g, mode = "out") / (group_size - 1)
    
    g_undir <- igraph::as_undirected(g, mode = "collapse", edge.attr.comb = "sum")
    ec <- tryCatch(
      igraph::eigen_centrality(g_undir)$vector,
      error = function(e) rep(NA, n_nodes)
    )
    eigen_mat[k, ] <- ec
  }
  
  list(
    in_strength  = in_strength_mat,
    out_strength = out_strength_mat,
    in_degree    = in_degree_mat,
    out_degree   = out_degree_mat,
    eigenvector  = eigen_mat
  )
}


# ══════════════════════════════════════════════════════════════════════
# Main function: extract grooming network metrics for immigrant males
# ══════════════════════════════════════════════════════════════════════
extract_male_grooming_metrics <- function(grooming_data, affiliative,
                                          emigrants, group_presence_matrices,
                                          num_draws = 1000,
                                          prior_a = 0.1, prior_b = 0.1) {
  
  # Prepare data
  grooming_data$date <- as.Date(grooming_data$date)
  grooming_data$idindividual1 <- tolower(as.character(grooming_data$idindividual1))
  grooming_data$idindividual2 <- tolower(as.character(grooming_data$idindividual2))
  
  affiliative$date <- as.Date(affiliative$date)
  
  results <- list()
  posterior_draws <- list()
  
  for (i in seq_len(nrow(emigrants))) {
    
    male_code  <- tolower(emigrants$Code[i])
    group_name <- emigrants$ImmigrationGp1[i]
    if (!is.na(emigrants$DateImmigration1[i])) {
      start_date <- as.Date(emigrants$DateImmigration1[i])
    } else if (!is.na(emigrants$DepartureNatalGp[i])) {
      start_date <- as.Date(emigrants$DepartureNatalGp[i])
    } else {
      start_date <- as.Date(NA)
    }
    end_date <- if (!is.na(start_date)) start_date + 365 else as.Date(NA)
    
    if (any(is.na(c(male_code, group_name, start_date)))) {
      message("Skipping ", male_code, " - missing required data")
      next
    }
    
    if (!group_name %in% names(group_presence_matrices)) {
      message("No presence matrix for group: ", group_name, " (", male_code, ")")
      next
    }
    
    # ── 1. Get observation days (from affiliative data) ────────────────────
    affiliative_period <- affiliative %>%
      filter(group_Abbreviation == group_name,
             date >= start_date,
             date <= end_date)
    
    days_observed <- length(unique(affiliative_period$date))
    
    if (days_observed == 0) {
      message("No observation days for ", male_code, " in ", group_name)
      next
    }
    
    # ── 2. Filter grooming data to group and time window ──────────────────
    group_groom <- grooming_data %>%
      filter(group_Abbreviation == group_name,
             date >= start_date,
             date <= end_date,
             !is.na(idindividual1), idindividual1 != "",
             !is.na(idindividual2), idindividual2 != "")
    
    if (nrow(group_groom) == 0) {
      message("No grooming data for ", male_code, " in ", group_name)
      next
    }
    
    # ── 3. Build directed dyadic counts ───────────────────────────────────
    # event = number of distinct days groomer groomed groomee
    # duration = total observation days in that period
    dyadic_data <- group_groom %>%
      group_by(groomer = idindividual1, groomee = idindividual2) %>%
      summarise(event = n_distinct(date), .groups = "drop") %>%
      filter(groomer != groomee) %>%
      mutate(duration = days_observed)
    
    if (nrow(dyadic_data) == 0) {
      message("No grooming dyads for ", male_code)
      next
    }
    
    all_nodes <- unique(c(dyadic_data$groomer, dyadic_data$groomee))
    
    if (!male_code %in% all_nodes) {
      message("Male ", male_code, " not in grooming network")
      next
    }
    
    # ── 4. Get median group size from presence matrix ─────────────────────
    presence <- group_presence_matrices[[group_name]]
    presence$Date <- as.Date(presence$Date)
    presence_period <- presence %>%
      filter(Date >= start_date, Date <= end_date)
    
    median_group_size <- NA
    if (nrow(presence_period) > 0) {
      daily_sizes <- rowSums(presence_period[, -1, drop = FALSE], na.rm = TRUE)
      median_group_size <- median(daily_sizes, na.rm = TRUE)
    }
    
    # ── 5. Compute Bayesian metrics ───────────────────────────────────────
    tryCatch({
      
      metrics <- compute_directed_metrics_from_posterior(
        dyadic_data = dyadic_data,
        all_nodes   = all_nodes,
        num_draws   = num_draws,
        prior_a     = prior_a,
        prior_b     = prior_b,
        group_size  = median_group_size
      )
      
      mi <- which(colnames(metrics$in_strength) == male_code)
      
      # Save full posterior draws
      posterior_draws[[male_code]] <- list(
        groom_in_strength  = metrics$in_strength[, mi],
        groom_out_strength = metrics$out_strength[, mi],
        groom_in_degree    = metrics$in_degree[, mi],
        groom_out_degree   = metrics$out_degree[, mi]
      )
      
      if (length(mi) == 0) {
        message("Male ", male_code, " not found in metric output")
        next
      }
      
      # Classic grooming rate (non-Bayesian)
      focal_groom <- group_groom %>%
        filter(idindividual1 == male_code | idindividual2 == male_code)
      grooming_rate <- length(unique(focal_groom$date)) / days_observed
      
      results[[i]] <- tibble(
        Code                       = male_code,
        group                      = group_name,
        start_date                 = start_date,
        end_date                   = end_date,
        # Bayesian in-strength (being groomed)
        in_strength_mean           = mean(metrics$in_strength[, mi]),
        in_strength_lower          = quantile(metrics$in_strength[, mi], 0.025),
        in_strength_upper          = quantile(metrics$in_strength[, mi], 0.975),
        # Bayesian out-strength (grooming others)
        out_strength_mean          = mean(metrics$out_strength[, mi]),
        out_strength_lower         = quantile(metrics$out_strength[, mi], 0.025),
        out_strength_upper         = quantile(metrics$out_strength[, mi], 0.975),
        # Bayesian in-degree
        in_degree_mean             = mean(metrics$in_degree[, mi]),
        in_degree_lower            = quantile(metrics$in_degree[, mi], 0.025),
        in_degree_upper            = quantile(metrics$in_degree[, mi], 0.975),
        # Bayesian out-degree
        out_degree_mean            = mean(metrics$out_degree[, mi]),
        out_degree_lower           = quantile(metrics$out_degree[, mi], 0.025),
        out_degree_upper           = quantile(metrics$out_degree[, mi], 0.975),
        # Bayesian eigenvector centrality (on undirected collapse)
        eigenvector_mean           = mean(metrics$eigenvector[, mi], na.rm = TRUE),
        eigenvector_lower          = quantile(metrics$eigenvector[, mi], 0.025, na.rm = TRUE),
        eigenvector_upper          = quantile(metrics$eigenvector[, mi], 0.975, na.rm = TRUE),
        # Classic metrics
        grooming_rate              = grooming_rate,
        # Sampling info
        days_observed              = days_observed,
        median_group_size          = median_group_size,
        n_dyads                    = nrow(dyadic_data)
      )
      
      message("Done: ", male_code, " (", group_name, ")")
      
    }, error = function(e) {
      message("Error for male ", male_code, ": ", e$message)
      NULL
    })
  }
  
  summaries <- bind_rows(results)
  list(
    summaries  = summaries,
    posteriors = posterior_draws
  )
}

groom_results <- extract_male_grooming_metrics(
  grooming_data, affiliative, emigrants, group_presence_matrices,
  num_draws = 1000
)

write.csv(groom_results$summaries,
          "/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/Pooja_chapter_matrank/Data/network_mets_groom_bayes.csv",
          row.names = FALSE)

saveRDS(groom_results$posteriors,
        "/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/Pooja_chapter_matrank/Data/groom_posteriors.rds")

cat("Grooming done.", length(groom_results$posteriors), "males with posteriors saved.\n")