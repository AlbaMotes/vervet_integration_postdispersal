#Script to calcualte sex ratios at the time of immigration
#Input files
#emigrants<-read.csv("/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Pooja_chapter_matrank/Data/All_potential_immigrants.csv")
# load life history IVP, conflict data and group_presence matrices
#load("/Users/alba/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_UNIL/Pooja_chapter_matrank/Data/Immigrants_environment.RData")

calculate_sex_ratio_at_immigration <- function(emigrants, life_hist, group_presence_matrices) {
  
  # Prepare life_hist data - get sex information
  life_hist_sex <- life_hist %>%
    select(Code, Sex) %>%
    mutate(
      Code_lower = tolower(as.character(Code)),
      Sex = toupper(as.character(Sex))
    )
  
  # Initialize result columns
  emigrants$total_group_size <- NA
  emigrants$sex_ratio_mf <- NA
  emigrants$prop_male <- NA
  emigrants$prop_unknown_sex <- NA
  
  # Process each emigrant
  for (i in 1:nrow(emigrants)) {
    
    emigrant_code <- emigrants$Code[i]
    group_name <- emigrants$ImmigrationGp1[i]
    immigration_date <- as.Date(emigrants$DateImmigration1[i])
    
    # Skip if missing data
    if (is.na(group_name) || is.na(immigration_date)) {
      next
    }
    
    # Check if group exists in presence matrices
    if (!group_name %in% names(group_presence_matrices)) {
      message(paste("Group", group_name, "not found in presence matrices for", emigrant_code))
      next
    }
    
    # Get presence matrix for this group
    presence_matrix <- group_presence_matrices[[group_name]]
    presence_matrix$Date <- as.Date(presence_matrix$Date)
    
    # Find the row closest to immigration date
    date_diffs <- abs(as.numeric(presence_matrix$Date - immigration_date))
    closest_date_idx <- which.min(date_diffs)
    
    # Get individuals present on that date
    presence_row <- presence_matrix[closest_date_idx, ]
    present_individuals <- names(presence_row)[presence_row == 1 & names(presence_row) != "Date"]
    
    if (length(present_individuals) == 0) {
      message(paste("No individuals present for", emigrant_code, "on", immigration_date))
      next
    }
    
    # Total group size
    total_size <- length(present_individuals)
    
    # Match with life history to get sex
    present_individuals_lower <- tolower(present_individuals)
    
    sex_info <- life_hist_sex %>%
      filter(Code_lower %in% present_individuals_lower)
    
    # Count by sex
    n_males <- sum(sex_info$Sex == "M", na.rm = TRUE)
    n_females <- sum(sex_info$Sex == "F", na.rm = TRUE)
    n_known_sex <- nrow(sex_info)
    n_unknown <- total_size - n_known_sex
    
    # Store results
    emigrants$total_group_size[i] <- total_size
    emigrants$prop_male[i] <- n_males / total_size
    emigrants$prop_unknown_sex[i] <- n_unknown / total_size
    
    # Calculate sex ratio (males/females)
    if (n_females > 0) {
      emigrants$sex_ratio_mf[i] <- n_males / n_females
    } else {
      emigrants$sex_ratio_mf[i] <- NA
    }
  }
  
  return(emigrants)
}


