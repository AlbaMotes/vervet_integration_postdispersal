library(dplyr)
library(tidyr)
library(igraph)
library(tibble)

source("/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/Pooja_chapter_matrank/Analyses/scan_extraction_1mandarm.R")

emigrants <- read.csv("/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/Pooja_chapter_matrank/Data/All_potential_immigrants_Josie_version.csv")
scans <- read.csv("/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/IVP_data/scans/scans_for_bison.csv")


# ══════════════════════════════════════════════════════════════════════
# Helper: compute network metrics from Beta posterior samples
# Returns both summaries AND full posterior draws per focal male
# ══════════════════════════════════════════════════════════════════════
compute_metrics_from_posterior <- function(dyadic_data, all_nodes, focal_node,
                                           num_draws = 1000,
                                           prior_a = 0.1, prior_b = 0.1,
                                           edge_thresh = 0.0001,
                                           group_size = NA) {
  
  n_dyads <- nrow(dyadic_data)
  n_nodes <- length(all_nodes)
  
  # Sample edge weights from Beta posterior
  draw_mat <- matrix(NA, nrow = n_dyads, ncol = num_draws)
  for (d in 1:n_dyads) {
    draw_mat[d, ] <- rbeta(num_draws,
                           prior_a + dyadic_data$event[d],
                           prior_b + dyadic_data$duration[d] - dyadic_data$event[d])
  }
  
  node1 <- as.character(dyadic_data$node_1)
  node2 <- as.character(dyadic_data$node_2)
  
  # Vectors to store focal node's metrics per draw
  strength_draws <- numeric(num_draws)
  degree_draws   <- numeric(num_draws)
  
  for (k in 1:num_draws) {
    
    adj <- matrix(0, n_nodes, n_nodes, dimnames = list(all_nodes, all_nodes))
    for (r in 1:n_dyads) {
      w <- draw_mat[r, k]
      if (w < edge_thresh) w <- 0
      adj[node1[r], node2[r]] <- w
      adj[node2[r], node1[r]] <- w
    }
    
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
    
    strength_draws[k] <- igraph::strength(g)[focal_node]
    # Normalized the same way as grooming degree in bison_groom_net_script.R:
    # proportion of the group (not just network nodes) the focal male is
    # connected to, rather than a raw unbounded edge count.
    degree_draws[k]   <- igraph::degree(g)[focal_node] / (group_size - 1)
  }
  
  list(
    strength_draws = strength_draws,
    degree_draws   = degree_draws
  )
}


# ══════════════════════════════════════════════════════════════════════
# Main function — saves full posteriors for BISoN-style analysis
# ══════════════════════════════════════════════════════════════════════
extract_male_network_metrics <- function(scans, emigrants, group_presence_matrices, num_draws = 1000) {
  
  contact_behaviours <- c(
    "allogrooming: groom", "allogrooming: being groomed",
    "affiliative: sit in contact", "affiliative: mouth to mouth",
    "affiliative: nurse", "affiliative: infant handling",
    "groom", "being groomed", "mouth to mouth", "mouth fed from",
    "feed from mouth", "nurse", "infant handle", "play-mount",
    "gr", "bgr", "mc", "bmc", "pm", "bpm", "sm", "bsm",
    "to", "bto", "ca", "em", "nu", "ih", "ii", "sg",
    "bite", "being bitten", "hit", "being hit", "grab",
    "hand on head", "hand on top of head (hh)", "attack", "chase",
    "fm", "bfm", "gb", "bgb", "bi", "bbi", "hi", "bhi",
    "fi", "hh", "bhh", "so", "bso",
    "sexual: mount", "sexual: be-mounted", "sexual: being mounted",
    "sexual: inspect", "sexual: solicit", "sexual: be solicited",
    "play: sexual", "mu", "bmu", "hg", "bhg", "bs", "bbs"
  )
  
  results <- list()
  # Store full posterior draws: list of lists, keyed by male Code
  posterior_draws <- list()
  
  for (i in seq_len(nrow(emigrants))) {
    
    male_code  <- emigrants$Code[i]
    group_name <- emigrants$ImmigrationGp1[i]
    start_date <- if (!is.na(emigrants$DateImmigration1[i])) {
      as.Date(emigrants$DateImmigration1[i])
    } else {
      as.Date(emigrants$DepartureNatalGp[i])
    }
    end_date <- if (!is.na(start_date)) start_date + 365 else as.Date(NA)
    
    if (any(is.na(c(male_code, group_name, start_date)))) {
      message("Skipping ", male_code, " - missing required data")
      next
    }
    
    scans_filtered <- scans %>%
      mutate(date = as.Date(date)) %>%
      filter(group_Abbreviation == group_name,
             !is.na(date),
             date >= start_date,
             date <= end_date)
    
    if (nrow(scans_filtered) == 0) {
      message("No scan data for male ", male_code)
      next
    }
    
    scans_filtered$year <- as.numeric(format(scans_filtered$date, "%Y"))
    proximity_data <- extract_proximity_pairs(scans_filtered)
    
    if (is.null(proximity_data) || nrow(proximity_data) == 0) {
      message("No proximity data for male ", male_code)
      next
    }
    
    from_focal <- scans_filtered %>%
      filter(!is.na(idindividual1) & idindividual1 != "") %>%
      mutate(
        is_contact    = trimws(tolower(behaviour))     %in% contact_behaviours |
          trimws(tolower(behaviourtype)) %in% contact_behaviours,
        has_proximity = (!is.na(ind1m)        & ind1m        != "") |
          (!is.na(indarmlength) & indarmlength != "") |
          (!is.na(nnadult)      & nnadult      != "") |
          (!is.na(nnjuvenile)   & nnjuvenile   != ""),
        exclude       = is_contact & !has_proximity
      ) %>%
      filter(!exclude) %>%
      distinct(idindividual1, scan_id) %>%
      rename(individual = idindividual1)
    
    from_proximity <- proximity_data %>%
      filter(!is.na(interaction_partner) & interaction_partner != "") %>%
      distinct(interaction_partner, scan_id) %>%
      rename(individual = interaction_partner)
    
    individual_scan_counts <- bind_rows(from_focal, from_proximity) %>%
      distinct() %>%
      group_by(individual) %>%
      summarise(scans_sampled = n_distinct(scan_id), .groups = "drop")
    
    dyadic_data <- proximity_data %>%
      filter(!is.na(idindividual1) & idindividual1 != "",
             !is.na(interaction_partner) & interaction_partner != "") %>%
      mutate(
        node_1 = pmin(idindividual1, interaction_partner),
        node_2 = pmax(idindividual1, interaction_partner)
      ) %>%
      filter(node_1 != node_2) %>%
      group_by(node_1, node_2) %>%
      summarise(event = n_distinct(scan_id), .groups = "drop") %>%
      left_join(individual_scan_counts, by = c("node_1" = "individual")) %>%
      rename(scans_node_1 = scans_sampled) %>%
      left_join(individual_scan_counts, by = c("node_2" = "individual")) %>%
      rename(scans_node_2 = scans_sampled) %>%
      mutate(
        duration = scans_node_1 + scans_node_2 - event,
        sri      = event / duration
      )
    
    if (nrow(dyadic_data) == 0 ||
        !male_code %in% c(dyadic_data$node_1, dyadic_data$node_2)) {
      message("Male ", male_code, " has no dyads in network")
      next
    }
    
    all_nodes <- unique(c(dyadic_data$node_1, dyadic_data$node_2))
    
    # ── Get median group size from presence matrix ─────────────────────────
    # Identical to bison_groom_net_script.R's step 4, so degree is
    # normalized the same way in both scripts.
    if (!group_name %in% names(group_presence_matrices)) {
      message("No presence matrix for group: ", group_name, " (", male_code, ")")
      next
    }
    presence <- group_presence_matrices[[group_name]]
    presence$Date <- as.Date(presence$Date)
    presence_period <- presence %>%
      filter(Date >= start_date, Date <= end_date)
    
    median_group_size <- NA
    if (nrow(presence_period) > 0) {
      daily_sizes <- rowSums(presence_period[, -1, drop = FALSE], na.rm = TRUE)
      median_group_size <- median(daily_sizes, na.rm = TRUE)
    }
    
    if (is.na(median_group_size) || median_group_size <= 1) {
      message("No usable group size for male ", male_code, " - skipping")
      next
    }
    
    tryCatch({
      
      metrics <- compute_metrics_from_posterior(
        dyadic_data = dyadic_data,
        all_nodes   = all_nodes,
        focal_node  = male_code,
        num_draws   = num_draws,
        prior_a     = 0.1,
        prior_b     = 0.1,
        group_size  = median_group_size
      )
      
      # Save full posterior draws for this male
      posterior_draws[[male_code]] <- list(
        prox_strength = metrics$strength_draws,
        prox_degree   = metrics$degree_draws
      )
      
      # Classic SRI-based strength
      male_dyads     <- dyadic_data %>%
        filter(node_1 == male_code | node_2 == male_code)
      strength_assoc <- sum(male_dyads$sri)
      
      results[[i]] <- tibble(
        Code             = male_code,
        group            = group_name,
        start_date       = start_date,
        end_date         = end_date,
        strength_mean    = mean(metrics$strength_draws),
        strength_lower   = quantile(metrics$strength_draws, 0.025),
        strength_upper   = quantile(metrics$strength_draws, 0.975),
        degree_mean      = mean(metrics$degree_draws),
        degree_lower     = quantile(metrics$degree_draws, 0.025),
        degree_upper     = quantile(metrics$degree_draws, 0.975),
        strength_assoc   = strength_assoc,
        n_dyads          = nrow(dyadic_data),
        median_group_size = median_group_size,
        n_scans          = individual_scan_counts$scans_sampled[
          individual_scan_counts$individual == male_code]
      )
      
      message("Done: ", male_code, " (", group_name, ")")
      
    }, error = function(e) {
      message("Metric extraction error for male ", male_code, ": ", e$message)
      NULL
    })
  }
  
  summaries <- bind_rows(results)
  
  # Return both summaries and posteriors
  list(
    summaries  = summaries,
    posteriors = posterior_draws
  )
}

# ══════════════════════════════════════════════════════════════════════
# Run and save
# ══════════════════════════════════════════════════════════════════════
group_presence_matrices <- readRDS("/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/Pooja_chapter_matrank/Data/group_presence_matrices.rds")

prox_results <- extract_male_network_metrics(scans, emigrants, group_presence_matrices, num_draws = 1000)

write.csv(prox_results$summaries,
          "/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/Pooja_chapter_matrank/Data/network_mets_prox_bayes.csv",
          row.names = FALSE)

saveRDS(prox_results$posteriors,
        "/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Vervet_projects/Pooja_chapter_matrank/Data/prox_posteriors.rds")
cat("Proximity done.", length(prox_results$posteriors), "males with posteriors saved.\n")