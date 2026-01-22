#calculate median group sizes to use in degree normalization


calculate_median_group_size <- function(potential_immig, group_presence_matrices, end_date_col = "elo_rating_date_all") {
  # Validate end_date_col parameter - check if column exists in the dataframe
  if (!end_date_col %in% colnames(potential_immig)) {
    stop(paste("Column", end_date_col, "not found in potential_immig. Available columns:", 
               paste(colnames(potential_immig), collapse = ", ")))
  }
  
  # Initialize results dataframe
  group_size_results <- data.frame(
    immigrant_code = character(),
    group = character(),
    median_group_size = numeric(),
    min_group_size = numeric(),
    max_group_size = numeric(),
    days_observed = numeric(),
    DateImmigration1 = as.Date(character()),
    end_date = as.Date(character()),
    end_date_type = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each immigrant
  for (i in 1:nrow(potential_immig)) {
    immigrant_code <- potential_immig$Code[i]
    group_abbreviation <- potential_immig$ImmigrationGp1[i]
    start_date <- potential_immig$DateImmigration1[i]
    end_date <- potential_immig[[end_date_col]][i]
    
    # Skip if any required data is missing
    if (is.na(immigrant_code) || is.na(group_abbreviation) || is.na(start_date) || is.na(end_date)) {
      warning(paste("Skipping immigrant", immigrant_code, "due to missing data"))
      next
    }
    
    # Ensure dates are in Date format
    start_date <- as.Date(start_date)
    end_date <- as.Date(end_date)
    
    # Check if end_date is after start_date
    if (end_date <= start_date) {
      warning(paste("Skipping immigrant", immigrant_code, "- end date is not after DateImmigration1"))
      next
    }
    
    # Get the presence matrix for this group using the abbreviation
    presence_matrix <- group_presence_matrices[[group_abbreviation]]
    if (is.null(presence_matrix)) {
      warning(paste("No presence matrix found for group", group_abbreviation, "for immigrant", immigrant_code))
      cat("Available group matrices:", paste(names(group_presence_matrices), collapse = ", "), "\n")
      next
    }
    
    # Convert row names (dates) to Date format if they aren't already
    matrix_dates <- as.Date(presence_matrix$Date)
    
    # Find rows (dates) within the period
    period_rows <- which(matrix_dates >= start_date & matrix_dates <= end_date)
    
    if (length(period_rows) == 0) {
      warning(paste("No presence data found for the period for immigrant", immigrant_code))
      next
    }
    
    # Extract the subset of the presence matrix for the relevant period
    period_matrix <- presence_matrix[period_rows, -1, drop = FALSE]  # Exclude Date column
    
    # Calculate group size for each day (count individuals present = 1, excluding the immigrant if present)
    daily_group_sizes <- c()
    
    for (row_idx in 1:nrow(period_matrix)) {
      total_present <- sum(period_matrix[row_idx, ] == 1, na.rm = TRUE)
      
      # Check if the immigrant is in the matrix and was present on this day
      immigrant_code_lower <- tolower(immigrant_code)
      immigrant_cols <- which(tolower(colnames(period_matrix)) == immigrant_code_lower)
      
      if (length(immigrant_cols) > 0) {
        # If immigrant is in the matrix, subtract 1 if they were present
        immigrant_present <- sum(period_matrix[row_idx, immigrant_cols] == 1, na.rm = TRUE)
        group_size <- total_present - immigrant_present
      } else {
        # If immigrant is not in the matrix, subtract 1 (assuming they were present)
        group_size <- total_present - 1
      }
      
      # Ensure group size is not negative
      group_size <- max(0, group_size)
      daily_group_sizes <- c(daily_group_sizes, group_size)
    }
    
    # Calculate median group size
    median_group_size <- median(daily_group_sizes, na.rm = TRUE)
    min_group_size <- min(daily_group_sizes, na.rm = TRUE)
    max_group_size <- max(daily_group_sizes, na.rm = TRUE)
    
    # Store results
    group_size_results <- rbind(group_size_results, data.frame(
      immigrant_code = immigrant_code,
      group = group_abbreviation,
      median_group_size = median_group_size,
      min_group_size = min_group_size,
      max_group_size = max_group_size,
      days_observed = length(daily_group_sizes),
      DateImmigration1 = start_date,
      end_date = end_date,
      end_date_type = end_date_col,
      stringsAsFactors = FALSE
    ))
    
    cat("Processed immigrant", immigrant_code, "- Median group size:", round(median_group_size, 1), 
        "over", length(daily_group_sizes), "days (using", end_date_col, ")\n")
  }
  
  return(group_size_results)
}
