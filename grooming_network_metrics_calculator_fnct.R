calculate_groom_network_metrics <- function(potential_immig, affiliative, grooming_data, group_presence_matrices, end_date_col = "elo_rating_date_one_year_all") {
  
  # Validate end_date_col parameter
  if (!end_date_col %in% names(potential_immig)) {
    stop(paste("Column", end_date_col, "not found in potential_immig dataset"))
  }
  
  # Clean presence matrices once - remove columns with NA or empty names
  group_presence_matrices <- lapply(group_presence_matrices, function(mat) {
    valid_cols <- !is.na(colnames(mat)) & colnames(mat) != "" & colnames(mat) != " "
    mat <- mat[, valid_cols, drop = FALSE]
    mat$Date <- as.Date(mat$Date)  # Convert dates once
    return(mat)
  })
  
  # Prepare potential_immig data
  potential_immig$Code_lower <- tolower(as.character(potential_immig$Code))
  potential_immig$DateImmigration1 <- as.Date(potential_immig$DateImmigration1)
  potential_immig[[end_date_col]] <- as.Date(potential_immig[[end_date_col]])
  
  # Filter to valid individuals (have group in presence matrices and valid end date)
  valid_individuals <- potential_immig %>%
    filter(
      ImmigrationGp1 %in% names(group_presence_matrices),
      !is.na(.data[[end_date_col]])
    )
  
  cat("Processing", nrow(valid_individuals), "individuals with valid data out of", nrow(potential_immig), "\n")
  
  # Prepare grooming data once
  grooming_data$date <- as.Date(grooming_data$date)
  grooming_data$idindividual1 <- tolower(as.character(grooming_data$idindividual1))
  grooming_data$idindividual2 <- tolower(as.character(grooming_data$idindividual2))
  grooming_data <- grooming_data[!is.na(grooming_data$idindividual1) & 
                                   !is.na(grooming_data$idindividual2) &
                                   !is.na(grooming_data$date), ]
  
  # Prepare affiliative data once
  affiliative$date <- as.Date(affiliative$date)
  affiliative <- affiliative[!is.na(affiliative$date), ]
  
  # Initialize results list (faster than rbind in loop)
  results_list <- vector("list", nrow(valid_individuals))
  
  # Process each valid individual
  for (i in 1:nrow(valid_individuals)) {
    
    male_code <- valid_individuals$Code_lower[i]
    group_name <- valid_individuals$ImmigrationGp1[i]
    start_date <- valid_individuals$DateImmigration1[i]
    end_date <- valid_individuals[[end_date_col]][i]
    
    # Get presence matrix for this group
    presence_matrix <- group_presence_matrices[[group_name]]
    
    # Filter presence matrix for tenure period
    presence_period <- presence_matrix[
      presence_matrix$Date >= start_date & presence_matrix$Date <= end_date,
    ]
    
    if (nrow(presence_period) == 0) {
      message(paste("Individual", male_code, "has no presence data in tenure period"))
      next
    }
    
    # Calculate median group size (vectorized)
    presence_data <- presence_period[, -1, drop = FALSE]  # Exclude Date column
    daily_group_sizes <- rowSums(presence_data, na.rm = TRUE)
    median_group_size <- median(daily_group_sizes, na.rm = TRUE)
    
    # Filter affiliative data for this group and time period
    affiliative_period <- affiliative[
      affiliative$group_Abbreviation == group_name &
        affiliative$date >= start_date &
        affiliative$date <= end_date,
    ]
    days_observed <- length(unique(affiliative_period$date))
    
    # Filter grooming data for this group and time period
    group_data <- grooming_data[
      grooming_data$group_Abbreviation == group_name &
        grooming_data$date >= start_date &
        grooming_data$date <= end_date,
    ]
    
    # Calculate grooming metrics (vectorized comparisons)
    focal_mask <- group_data$idindividual1 == male_code | group_data$idindividual2 == male_code
    in_mask <- group_data$idindividual2 == male_code
    out_mask <- group_data$idindividual1 == male_code
    
    total_grooming_days <- length(unique(group_data$date[focal_mask]))
    in_grooming_days <- length(unique(group_data$date[in_mask]))
    out_grooming_days <- length(unique(group_data$date[out_mask]))
    
    grooming_rate <- if(days_observed > 0) total_grooming_days / days_observed else 0
    in_grooming_rate <- if(days_observed > 0) in_grooming_days / days_observed else 0
    out_grooming_rate <- if(days_observed > 0) out_grooming_days / days_observed else 0
    
    # Network metrics
    in_degree_normalized <- 0
    out_degree_normalized <- 0
    total_degree_normalized <- 0
    eigenvector_centrality <- 0
    
    if (nrow(group_data) > 0) {
      
      # Create edge list with unique date counts per dyad
      edge_list <- group_data %>%
        group_by(idindividual1, idindividual2) %>%
        summarise(weight = n_distinct(date), .groups = 'drop')
      
      # Create directed graph
      g <- graph_from_data_frame(edge_list, directed = TRUE)
      
      # Check if focal individual is in the network
      if (male_code %in% V(g)$name) {
        
        # Calculate degrees
        in_deg <- igraph::degree(g, male_code, mode = "in")
        out_deg <- igraph::degree(g, male_code, mode = "out")
        
        # Count unique partners
        unique_partners <- length(unique(c(
          igraph::neighbors(g, male_code, mode = "in"),
          igraph::neighbors(g, male_code, mode = "out")
        )))
        
        # Normalize by available partners
        available_partners <- median_group_size - 1
        if (available_partners > 0) {
          in_degree_normalized <- in_deg / available_partners
          out_degree_normalized <- out_deg / available_partners
          total_degree_normalized <- unique_partners / available_partners
        }
        
        # Calculate eigenvector centrality
        tryCatch({
          g_undirected <- as.undirected(g, mode = "collapse", edge.attr.comb = "sum")
          eigen_cent <- eigen_centrality(g_undirected, directed = FALSE, weights = E(g_undirected)$weight)
          
          if (male_code %in% V(g_undirected)$name) {
            eigenvector_centrality <- eigen_cent$vector[male_code]
          }
        }, error = function(e) {
          message(paste("Error calculating eigenvector centrality for", male_code, ":", e$message))
        })
      }
    }
    
    # Store results in list
    results_list[[i]] <- data.frame(
      immigrant_code = male_code,
      group = group_name,
      grooming_rate = grooming_rate,
      in_grooming_rate = in_grooming_rate,
      out_grooming_rate = out_grooming_rate,
      in_degree_normalized = in_degree_normalized,
      out_degree_normalized = out_degree_normalized,
      total_degree_normalized = total_degree_normalized,
      eigenvector_centrality = eigenvector_centrality,
      total_grooming_days = total_grooming_days,
      in_grooming_days = in_grooming_days,
      out_grooming_days = out_grooming_days,
      days_observed = days_observed,
      median_group_size = median_group_size,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all results at once (much faster than rbind in loop)
  results <- bind_rows(results_list)
  
  # Add back individuals who were filtered out with NA values
  all_individuals <- data.frame(
    immigrant_code = potential_immig$Code_lower,
    stringsAsFactors = FALSE
  )
  
  results <- all_individuals %>%
    left_join(results, by = "immigrant_code")
  
  return(results)
}
