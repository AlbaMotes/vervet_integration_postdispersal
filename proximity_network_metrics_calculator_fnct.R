# Function to calculate network metrics for males using scan_list data 
#(output of extract_scans_fnct.R)

calculate_proximity_network_metrics <- function(scan_list, potential_immig) {
  # Initialize results dataframe
  network_metrics <- data.frame(
    immigrant_code = character(),
    group = character(),
    degree = numeric(),
    strength = numeric(),
    strength_assoc = numeric(),
    eigenvector_centrality = numeric(),
    total_scans = numeric(),
    total_dyads = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each immigrant male
  for (immigrant_code in names(scan_list)) {
    # Get the scan data for this immigrant (already processed with interaction_partner column)
    scan_data <- scan_list[[immigrant_code]]
    
    if (is.null(scan_data) || nrow(scan_data) == 0) {
      message(paste("No scan data for immigrant", immigrant_code, "- storing with NA metrics"))
      immigrant_info <- potential_immig[potential_immig$Code == immigrant_code, ]
      group_name <- if(nrow(immigrant_info) > 0) immigrant_info$ImmigrationGp1[1] else NA
      
      network_metrics <- rbind(network_metrics, data.frame(
        immigrant_code = immigrant_code,
        group = group_name,
        degree = NA,
        strength = NA,
        strength_assoc = NA,
        eigenvector_centrality = NA,
        total_scans = NA,
        total_dyads = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # Get immigrant info
    immigrant_info <- potential_immig[potential_immig$Code == immigrant_code, ]
    if (nrow(immigrant_info) == 0) {
      warning(paste("No immigrant info found for", immigrant_code))
      next
    }
    
    group_name <- immigrant_info$ImmigrationGp1[1]
    
    # Calculate individual scan counts for ALL individuals (both idindividual1 and interaction_partner)
    individual_scan_counts <- scan_data %>%
      filter(!is.na(idindividual1) & idindividual1 != "") %>%
      pivot_longer(cols = c(idindividual1, interaction_partner), 
                   values_to = "individual", 
                   names_to = "role") %>%
      filter(!is.na(individual) & individual != "") %>%
      group_by(individual) %>%
      summarise(scans_sampled = n_distinct(scan_id), .groups = 'drop')
    
    # Create standardized dyads (always put alphabetically smaller ID first)
    scan_dyads <- scan_data %>%
      filter(!is.na(interaction_partner) & interaction_partner != "" &
               !is.na(idindividual1) & idindividual1 != "") %>%
      mutate(
        individual_a = pmin(idindividual1, interaction_partner),
        individual_b = pmax(idindividual1, interaction_partner)
      ) %>%
      filter(individual_a != individual_b)  # Remove self-loops
    
    cat("Valid dyads for immigrant", immigrant_code, ":", nrow(scan_dyads), "\n")
    
    if (nrow(scan_dyads) == 0) {
      warning(paste("No valid dyads found for immigrant", immigrant_code))
      next
    }
    
    # Create edgelist by counting co-occurrences per scan (following past method)
    edgelist <- scan_dyads %>%
      group_by(individual_a, individual_b) %>%
      summarise(
        co_occurrence_scans = n_distinct(scan_id),  # Number of scans they appeared together
        total_observations = n(),                    # Total observations (for reference)
        .groups = 'drop'
      ) %>%
      # Add individual scan counts
      left_join(
        individual_scan_counts %>% select(individual, scans_individual_a = scans_sampled),
        by = c("individual_a" = "individual")
      ) %>%
      left_join(
        individual_scan_counts %>% select(individual, scans_individual_b = scans_sampled),
        by = c("individual_b" = "individual")
      ) %>%
      # Calculate association metrics (same as past method)
      mutate(
        total_dyad_scans = scans_individual_a + scans_individual_b,
        association_rate = co_occurrence_scans / total_dyad_scans
      ) %>%
      # Remove dyads where one individual has no scan count
      filter(!is.na(scans_individual_a) & !is.na(scans_individual_b))
    
    if (nrow(edgelist) == 0) {
      warning(paste("No valid edgelist created for immigrant", immigrant_code))
      next
    }
    
    # Create igraph object with both edge weights
    g <- graph_from_data_frame(
      d = edgelist[, c("individual_a", "individual_b", "co_occurrence_scans", "association_rate")],
      directed = FALSE
    )
    
    # Check if immigrant is in the network
    if (!immigrant_code %in% V(g)$name) {
      warning(paste("Immigrant", immigrant_code, "not found in network - storing with NA metrics"))
      
      network_metrics <- rbind(network_metrics, data.frame(
        immigrant_code = immigrant_code,
        group = group_name,
        degree = NA,
        strength = NA,
        strength_assoc = NA,
        eigenvector_centrality = NA,
        total_scans = NA,
        total_dyads = nrow(edgelist),
        stringsAsFactors = FALSE
      ))
      next
    }    
    # Calculate network metrics for the immigrant
    degree_val <- igraph::degree(g,  immigrant_code)
    strength_val <- strength(g, v = immigrant_code, weights = E(g)$co_occurrence_scans)  # Pure strength
    strength_assoc_val <- strength(g, v = immigrant_code, weights = E(g)$association_rate)  # Association rate strength
    
    # Calculate eigenvector centrality (handle potential errors)
    eigen_cent <- tryCatch({
      eigen_centrality(g, weights = E(g)$association_rate)$vector[immigrant_code]
    }, error = function(e) {
      warning(paste("Could not calculate eigenvector centrality for", immigrant_code, ":", e$message))
      NA
    })
    
    # Get total scans for the immigrant
    immigrant_scans <- individual_scan_counts$scans_sampled[individual_scan_counts$individual == immigrant_code]
    if (length(immigrant_scans) == 0) immigrant_scans <- NA
    
    # Store results
    network_metrics <- rbind(network_metrics, data.frame(
      immigrant_code = immigrant_code,
      group = group_name,
      degree = degree_val,
      strength = strength_val,              # Sum of co_occurrence_scans
      strength_assoc = strength_assoc_val,  # Sum of association rates
      eigenvector_centrality = ifelse(is.null(eigen_cent), NA, eigen_cent),
      total_scans = immigrant_scans,
      total_dyads = nrow(edgelist),
      stringsAsFactors = FALSE
    ))
  }
  
  return(network_metrics)
}
