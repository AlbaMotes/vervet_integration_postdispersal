##Function to extract scan data for the desired time period for specific individuals

#The function uses as default for the end date of the metrics calculation the date the 
#elo rating is calculated. But it can be adapted to include all data collected 
#in the new group by specifying end_date_col = "LastSeen1"


extract_male_scan_interactions <- function(scans, potential_immig, end_date_col = "elo_rating_date_male_male") {
  # Create an empty list to store results for each male
  scan_list <- list()
  
  # Validate end_date_col parameter
  if (!end_date_col %in% names(potential_immig)) {
    stop(paste("Column", end_date_col, "not found in potential_immig dataset"))
  }
  
  # Process each male in the potential_immig dataset
  for (i in 1:nrow(potential_immig)) {
    male_code <- potential_immig$Code[i]
    scan_list[[male_code]] <- data.frame()
    group_name <- potential_immig$ImmigrationGp1[i]
    start_date <- potential_immig$DateImmigration1[i]
    end_date <- potential_immig[[end_date_col]][i]  # Use specified column
    
    # Skip if any required data is missing
    if (is.na(male_code) || is.na(group_name) || is.na(start_date) || is.na(end_date)) {
      next
    }
    
    # Ensure dates are in Date format
    start_date <- as.Date(start_date)
    end_date <- as.Date(end_date)
    
    # Convert scan date column to Date format if it's character
    scans$date <- as.Date(scans$date)
    
    # Filter scan data for the specific group and time period
    group_data <- scans[
      scans$group_Abbreviation == group_name &
        !is.na(scans$date) &
        scans$date >= start_date &
        scans$date <= end_date,
    ]
    
    # Skip if no data found
    if (nrow(group_data) == 0) {
      message(paste("Male", male_code, "has no scan data - keeping with empty data"))
      next
    }
    
    # Add year column for temporal logic
    group_data$year <- as.numeric(format(group_data$date, "%Y"))
    
    # Apply the correct hierarchy for proximity data
    proximity_data <- list()
    
    # Split data by year
    data_before_2022 <- group_data %>% filter(year < 2022)
    data_from_2022 <- group_data %>% filter(year >= 2022)
    
    # BEFORE 2022: 1m priority hierarchy
    if (nrow(data_before_2022) > 0) {
      processed_scans <- c()
      
      # 1. Check ind1m first
      if (any(!is.na(data_before_2022$ind1m) & data_before_2022$ind1m != "")) {
        ind1m_data <- data_before_2022 %>%
          filter(!is.na(ind1m) & ind1m != "") %>%
          separate_rows(ind1m, sep = ";") %>%
          mutate(interaction_partner = trimws(ind1m)) %>%
          filter(interaction_partner != "") %>%
          select(-ind1m)
        proximity_data <- append(proximity_data, list(ind1m_data))
        processed_scans <- unique(ind1m_data$scan_id)
      }
      
      # 2. Check distance-based 1m for remaining scans
      remaining_data <- data_before_2022 %>% filter(!scan_id %in% processed_scans)
      if (nrow(remaining_data) > 0) {
        distance_1m_data <- list()
        
        # nnadult when distancenna = 1m
        nnadult_1m <- remaining_data %>%
          filter(!is.na(distancenna) & distancenna == "1m" &
                   !is.na(nnadult) & nnadult != "") %>%
          separate_rows(nnadult, sep = ";") %>%
          mutate(interaction_partner = trimws(tolower(nnadult))) %>%
          filter(interaction_partner != "") %>%
          select(-nnadult)
        if (nrow(nnadult_1m) > 0) {
          distance_1m_data <- append(distance_1m_data, list(nnadult_1m))
        }
        
        # nnjuvenile when distancennj = 1m
        nnjuvenile_1m <- remaining_data %>%
          filter(!is.na(distancennj) & distancennj == "1m" &
                   !is.na(nnjuvenile) & nnjuvenile != "") %>%
          separate_rows(nnjuvenile, sep = ";") %>%
          mutate(interaction_partner = trimws(tolower(nnjuvenile))) %>%
          filter(interaction_partner != "") %>%
          select(-nnjuvenile)
        if (nrow(nnjuvenile_1m) > 0) {
          distance_1m_data <- append(distance_1m_data, list(nnjuvenile_1m))
        }
        
        if (length(distance_1m_data) > 0) {
          combined_1m <- bind_rows(distance_1m_data)
          proximity_data <- append(proximity_data, list(combined_1m))
          processed_scans <- c(processed_scans, unique(combined_1m$scan_id))
        }
      }
      
      # 3. Check indarmlength and Arm length distances for final remaining scans
      final_remaining <- data_before_2022 %>% filter(!scan_id %in% processed_scans)
      if (nrow(final_remaining) > 0) {
        arm_processed_scans <- c()
        
        # Check indarmlength
        if (any(!is.na(final_remaining$indarmlength) & final_remaining$indarmlength != "")) {
          indarmlength_data <- final_remaining %>%
            filter(!is.na(indarmlength) & indarmlength != "") %>%
            separate_rows(indarmlength, sep = ";") %>%
            mutate(interaction_partner = trimws(indarmlength)) %>%
            filter(interaction_partner != "") %>%
            select(-indarmlength)
          proximity_data <- append(proximity_data, list(indarmlength_data))
          arm_processed_scans <- unique(indarmlength_data$scan_id)
        }
        
        # Check Arm length distances for remaining scans
        arm_remaining <- final_remaining %>% filter(!scan_id %in% arm_processed_scans)
        if (nrow(arm_remaining) > 0) {
          distance_arm_data <- list()
          
          # nnadult when distancenna = Arm length
          nnadult_arm <- arm_remaining %>%
            filter(!is.na(distancenna) & distancenna == "Arm length" &
                     !is.na(nnadult) & nnadult != "") %>%
            separate_rows(nnadult, sep = ";") %>%
            mutate(interaction_partner = trimws(tolower(nnadult))) %>%
            filter(interaction_partner != "") %>%
            select(-nnadult)
          if (nrow(nnadult_arm) > 0) {
            distance_arm_data <- append(distance_arm_data, list(nnadult_arm))
          }
          
          # nnjuvenile when distancennj = Arm length
          nnjuvenile_arm <- arm_remaining %>%
            filter(!is.na(distancennj) & distancennj == "Arm length" &
                     !is.na(nnjuvenile) & nnjuvenile != "") %>%
            separate_rows(nnjuvenile, sep = ";") %>%
            mutate(interaction_partner = trimws(tolower(nnjuvenile))) %>%
            filter(interaction_partner != "") %>%
            select(-nnjuvenile)
          if (nrow(nnjuvenile_arm) > 0) {
            distance_arm_data <- append(distance_arm_data, list(nnjuvenile_arm))
          }
          
          if (length(distance_arm_data) > 0) {
            combined_arm <- bind_rows(distance_arm_data)
            proximity_data <- append(proximity_data, list(combined_arm))
          }
        }
      }
    }
    
    # FROM 2022 ONWARDS: Arm length priority hierarchy
    if (nrow(data_from_2022) > 0) {
      processed_scans_2022 <- c()
      
      # 1. Check indarmlength first
      if (any(!is.na(data_from_2022$indarmlength) & data_from_2022$indarmlength != "")) {
        indarmlength_data_2022 <- data_from_2022 %>%
          filter(!is.na(indarmlength) & indarmlength != "") %>%
          separate_rows(indarmlength, sep = ";") %>%
          mutate(interaction_partner = trimws(indarmlength)) %>%
          filter(interaction_partner != "") %>%
          select(-indarmlength)
        proximity_data <- append(proximity_data, list(indarmlength_data_2022))
        processed_scans_2022 <- unique(indarmlength_data_2022$scan_id)
      }
      
      # 2. Check distance-based Arm length for remaining scans
      remaining_data_2022 <- data_from_2022 %>% filter(!scan_id %in% processed_scans_2022)
      if (nrow(remaining_data_2022) > 0) {
        distance_arm_data_2022 <- list()
        
        # nnadult when distancenna = Arm length
        nnadult_arm_2022 <- remaining_data_2022 %>%
          filter(!is.na(distancenna) & distancenna == "Arm length" &
                   !is.na(nnadult) & nnadult != "") %>%
          separate_rows(nnadult, sep = ";") %>%
          mutate(interaction_partner = trimws(tolower(nnadult))) %>%
          filter(interaction_partner != "") %>%
          select(-nnadult)
        if (nrow(nnadult_arm_2022) > 0) {
          distance_arm_data_2022 <- append(distance_arm_data_2022, list(nnadult_arm_2022))
        }
        
        # nnjuvenile when distancennj = Arm length
        nnjuvenile_arm_2022 <- remaining_data_2022 %>%
          filter(!is.na(distancennj) & distancennj == "Arm length" &
                   !is.na(nnjuvenile) & nnjuvenile != "") %>%
          separate_rows(nnjuvenile, sep = ";") %>%
          mutate(interaction_partner = trimws(tolower(nnjuvenile))) %>%
          filter(interaction_partner != "") %>%
          select(-nnjuvenile)
        if (nrow(nnjuvenile_arm_2022) > 0) {
          distance_arm_data_2022 <- append(distance_arm_data_2022, list(nnjuvenile_arm_2022))
        }
        
        if (length(distance_arm_data_2022) > 0) {
          combined_arm_2022 <- bind_rows(distance_arm_data_2022)
          proximity_data <- append(proximity_data, list(combined_arm_2022))
          processed_scans_2022 <- c(processed_scans_2022, unique(combined_arm_2022$scan_id))
        }
      }
      
      # 3. Check ind1m and 1m distances for final remaining scans
      final_remaining_2022 <- data_from_2022 %>% filter(!scan_id %in% processed_scans_2022)
      if (nrow(final_remaining_2022) > 0) {
        one_m_processed_scans <- c()
        
        # Check ind1m
        if (any(!is.na(final_remaining_2022$ind1m) & final_remaining_2022$ind1m != "")) {
          ind1m_data_2022 <- final_remaining_2022 %>%
            filter(!is.na(ind1m) & ind1m != "") %>%
            separate_rows(ind1m, sep = ";") %>%
            mutate(interaction_partner = trimws(ind1m)) %>%
            filter(interaction_partner != "") %>%
            select(-ind1m)
          proximity_data <- append(proximity_data, list(ind1m_data_2022))
          one_m_processed_scans <- unique(ind1m_data_2022$scan_id)
        }
        
        # Check 1m distances for remaining scans
        one_m_remaining <- final_remaining_2022 %>% filter(!scan_id %in% one_m_processed_scans)
        if (nrow(one_m_remaining) > 0) {
          distance_1m_data_2022 <- list()
          
          # nnadult when distancenna = 1m
          nnadult_1m_2022 <- one_m_remaining %>%
            filter(!is.na(distancenna) & distancenna == "1m" &
                     !is.na(nnadult) & nnadult != "") %>%
            separate_rows(nnadult, sep = ";") %>%
            mutate(interaction_partner = trimws(tolower(nnadult))) %>%
            filter(interaction_partner != "") %>%
            select(-nnadult)
          if (nrow(nnadult_1m_2022) > 0) {
            distance_1m_data_2022 <- append(distance_1m_data_2022, list(nnadult_1m_2022))
          }
          
          # nnjuvenile when distancennj = 1m
          nnjuvenile_1m_2022 <- one_m_remaining %>%
            filter(!is.na(distancennj) & distancennj == "1m" &
                     !is.na(nnjuvenile) & nnjuvenile != "") %>%
            separate_rows(nnjuvenile, sep = ";") %>%
            mutate(interaction_partner = trimws(tolower(nnjuvenile))) %>%
            filter(interaction_partner != "") %>%
            select(-nnjuvenile)
          if (nrow(nnjuvenile_1m_2022) > 0) {
            distance_1m_data_2022 <- append(distance_1m_data_2022, list(nnjuvenile_1m_2022))
          }
          
          if (length(distance_1m_data_2022) > 0) {
            combined_1m_2022 <- bind_rows(distance_1m_data_2022)
            proximity_data <- append(proximity_data, list(combined_1m_2022))
          }
        }
      }
    }
    
    # Combine all proximity data
    if (length(proximity_data) > 0) {
      group_data_final <- bind_rows(proximity_data) %>%
        distinct()  # Remove any duplicates
      scan_list[[male_code]] <- group_data_final
    } else {
      message(paste("Male", male_code, "has no proximity data - keeping with empty data"))
      # scan_list[[male_code]] already initialized as empty dataframe at the start
    }
  }
  
  return(scan_list)
}
