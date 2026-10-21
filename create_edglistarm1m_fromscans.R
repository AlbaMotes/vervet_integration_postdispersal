library(tidyr)
library(dplyr)

create_scan_edgelist <- function(scan_filtered) {
  
  # Extract year from date
  scan_filtered$year <- as.numeric(format(as.Date(scan_filtered$date), "%Y"))
  
  # Apply the hierarchy for proximity data
  proximity_data <- list()
  
  # Split data by year
  data_before_2022 <- scan_filtered %>% filter(year < 2022)
  data_from_2022   <- scan_filtered %>% filter(year >= 2022)
  
  # BEFORE 2022: 1m priority hierarchy
  if (nrow(data_before_2022) > 0) {
    processed_scans <- c()
    
    # 1. Check ind1m first
    if (any(!is.na(data_before_2022$ind1m) & data_before_2022$ind1m != "")) {
      ind1m_data <- data_before_2022 %>%
        filter(!is.na(ind1m) & ind1m != "") %>%
        separate_rows(ind1m, sep = ";") %>%
        mutate(proximity_partner = trimws(ind1m)) %>%
        filter(proximity_partner != "") %>%
        dplyr::select(-ind1m)
      proximity_data  <- append(proximity_data, list(ind1m_data))
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
        mutate(proximity_partner = trimws(tolower(nnadult))) %>%
        filter(proximity_partner != "") %>%
        dplyr::select(-nnadult)
      if (nrow(nnadult_1m) > 0) {
        distance_1m_data <- append(distance_1m_data, list(nnadult_1m))
      }
      
      # nnjuvenile when distancennj = 1m
      nnjuvenile_1m <- remaining_data %>%
        filter(!is.na(distancennj) & distancennj == "1m" &
                 !is.na(nnjuvenile) & nnjuvenile != "") %>%
        separate_rows(nnjuvenile, sep = ";") %>%
        mutate(proximity_partner = trimws(tolower(nnjuvenile))) %>%
        filter(proximity_partner != "") %>%
        dplyr::select(-nnjuvenile)
      if (nrow(nnjuvenile_1m) > 0) {
        distance_1m_data <- append(distance_1m_data, list(nnjuvenile_1m))
      }
      
      if (length(distance_1m_data) > 0) {
        combined_1m     <- bind_rows(distance_1m_data)
        proximity_data  <- append(proximity_data, list(combined_1m))
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
          mutate(proximity_partner = trimws(indarmlength)) %>%
          filter(proximity_partner != "") %>%
          dplyr::select(-indarmlength)
        proximity_data      <- append(proximity_data, list(indarmlength_data))
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
          mutate(proximity_partner = trimws(tolower(nnadult))) %>%
          filter(proximity_partner != "") %>%
          dplyr::select(-nnadult)
        if (nrow(nnadult_arm) > 0) {
          distance_arm_data <- append(distance_arm_data, list(nnadult_arm))
        }
        
        # nnjuvenile when distancennj = Arm length
        nnjuvenile_arm <- arm_remaining %>%
          filter(!is.na(distancennj) & distancennj == "Arm length" &
                   !is.na(nnjuvenile) & nnjuvenile != "") %>%
          separate_rows(nnjuvenile, sep = ";") %>%
          mutate(proximity_partner = trimws(tolower(nnjuvenile))) %>%
          filter(proximity_partner != "") %>%
          dplyr::select(-nnjuvenile)
        if (nrow(nnjuvenile_arm) > 0) {
          distance_arm_data <- append(distance_arm_data, list(nnjuvenile_arm))
        }
        
        if (length(distance_arm_data) > 0) {
          combined_arm   <- bind_rows(distance_arm_data)
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
        mutate(proximity_partner = trimws(indarmlength)) %>%
        filter(proximity_partner != "") %>%
        dplyr::select(-indarmlength)
      proximity_data       <- append(proximity_data, list(indarmlength_data_2022))
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
        mutate(proximity_partner = trimws(tolower(nnadult))) %>%
        filter(proximity_partner != "") %>%
        dplyr::select(-nnadult)
      if (nrow(nnadult_arm_2022) > 0) {
        distance_arm_data_2022 <- append(distance_arm_data_2022, list(nnadult_arm_2022))
      }
      
      # nnjuvenile when distancennj = Arm length
      nnjuvenile_arm_2022 <- remaining_data_2022 %>%
        filter(!is.na(distancennj) & distancennj == "Arm length" &
                 !is.na(nnjuvenile) & nnjuvenile != "") %>%
        separate_rows(nnjuvenile, sep = ";") %>%
        mutate(proximity_partner = trimws(tolower(nnjuvenile))) %>%
        filter(proximity_partner != "") %>%
        dplyr::select(-nnjuvenile)
      if (nrow(nnjuvenile_arm_2022) > 0) {
        distance_arm_data_2022 <- append(distance_arm_data_2022, list(nnjuvenile_arm_2022))
      }
      
      if (length(distance_arm_data_2022) > 0) {
        combined_arm_2022    <- bind_rows(distance_arm_data_2022)
        proximity_data       <- append(proximity_data, list(combined_arm_2022))
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
          mutate(proximity_partner = trimws(ind1m)) %>%
          filter(proximity_partner != "") %>%
          dplyr::select(-ind1m)
        proximity_data        <- append(proximity_data, list(ind1m_data_2022))
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
          mutate(proximity_partner = trimws(tolower(nnadult))) %>%
          filter(proximity_partner != "") %>%
          dplyr::select(-nnadult)
        if (nrow(nnadult_1m_2022) > 0) {
          distance_1m_data_2022 <- append(distance_1m_data_2022, list(nnadult_1m_2022))
        }
        
        # nnjuvenile when distancennj = 1m
        nnjuvenile_1m_2022 <- one_m_remaining %>%
          filter(!is.na(distancennj) & distancennj == "1m" &
                   !is.na(nnjuvenile) & nnjuvenile != "") %>%
          separate_rows(nnjuvenile, sep = ";") %>%
          mutate(proximity_partner = trimws(tolower(nnjuvenile))) %>%
          filter(proximity_partner != "") %>%
          dplyr::select(-nnjuvenile)
        if (nrow(nnjuvenile_1m_2022) > 0) {
          distance_1m_data_2022 <- append(distance_1m_data_2022, list(nnjuvenile_1m_2022))
        }
        
        if (length(distance_1m_data_2022) > 0) {
          combined_1m_2022 <- bind_rows(distance_1m_data_2022)
          proximity_data   <- append(proximity_data, list(combined_1m_2022))
        }
      }
    }
  }
  
  # Combine all proximity data
  if (length(proximity_data) > 0) {
    scan_data_expanded <- bind_rows(proximity_data) %>%
      distinct()
  } else {
    warning("No proximity data found")
    return(data.frame())
  }
  
  # Exclude contact behaviours to count true 0s where individuals are present but not in proximity
  contact_behaviours <- c(
    # Affiliative - full labels
    "allogrooming: groom", "allogrooming: being groomed",
    "affiliative: sit in contact",
    "affiliative: mouth to mouth", "affiliative: nurse",
    "affiliative: infant handling",
    "groom", "being groomed",
    "mouth to mouth", "mouth fed from", "mouth fed from",
    "feed from mouth", "feed from mouth",
    "nurse",
    "infant handle",
    "play-mount",
    
    # Affiliative - abbreviations
    "gr", "bgr",
    "mc", "bmc",
    "pm", "bpm",
    "sm", "bsm",
    "to", "bto",
    "ca",
    "em",
    "nu",
    "ih", "ii",
    "sg",
    
    # Agonistic - full labels
    "bite", "being bitten",
    "hit", "being hit",
    "grab",
    "hand on head", "hand on top of head (hh)",
    "attack", "attack; be bitten",
    "chase",
    
    # Agonistic - abbreviations
    "fm", "bfm",
    "gb", "bgb",
    "bi", "bbi",
    "hi", "bhi",
    "fi",
    "hh", "bhh",
    "so", "bso",
    
    # Sexual - full labels
    "sexual: mount", "sexual: be-mounted", "sexual: being mounted",
    "sexual: inspect", "sexual: inspect",
    "sexual: solicit", "sexual: be solicited",
    "play: sexual",
    
    # Sexual - abbreviations
    "mu", "bmu",
    "hg", "bhg",
    "bs", "bbs"
  )
  
  # Calculate individual scan counts per group and season_id including true zeros were no individual was in proximity
  #but excluding cases where individuals were in physical contact.

  # Get all individuals from idindividual1 (with contact behaviour exclusion)
  from_focal <- scan_filtered %>%
    filter(!is.na(idindividual1) & idindividual1 != "") %>%
    mutate(
      is_contact_behaviour = trimws(tolower(behaviour))     %in% contact_behaviours |
        trimws(tolower(behaviourtype)) %in% contact_behaviours,
      has_proximity = (
        (!is.na(ind1m)        & ind1m        != "") |
          (!is.na(indarmlength) & indarmlength != "") |
          (!is.na(nnadult)      & nnadult      != "") |
          (!is.na(nnjuvenile)   & nnjuvenile   != "")
      ),
      exclude_scan = is_contact_behaviour & !has_proximity
    ) %>%
    filter(!exclude_scan) %>%
    distinct(idindividual1, group, season_id, scan_id) %>%
    rename(individual = idindividual1)
  
  # Get all individuals appearing as proximity partners
  # These by definition have proximity info recorded, so no exclusion needed
  from_proximity <- scan_data_expanded %>%
    filter(!is.na(proximity_partner) & proximity_partner != "") %>%
    distinct(proximity_partner, group, season_id, scan_id) %>%
    rename(individual = proximity_partner)
  
  # Combine and count
  individual_scan_counts <- bind_rows(from_focal, from_proximity) %>%
    distinct() %>%
    group_by(individual, group, season_id) %>%
    summarise(scans_sampled = n_distinct(scan_id), .groups = 'drop')

  
  # Create standardized dyads (always put alphabetically smaller ID first)
  scan_filtered_dyads <- scan_data_expanded %>%
    filter(!is.na(proximity_partner) & proximity_partner != "" &
             !is.na(idindividual1) & idindividual1 != "") %>%
    mutate(
      individual_a = pmin(idindividual1, proximity_partner),
      individual_b = pmax(idindividual1, proximity_partner)
    ) %>%
    filter(individual_a != individual_b)  # Remove self-loops
  
  # Create edgelist by counting co-occurrences per scan
  edgelist <- scan_filtered_dyads %>%
    group_by(individual_a, individual_b, group, season_id) %>%
    summarise(
      co_occurrence_scans = n_distinct(scan_id),
      total_observations  = n(),
      .groups = 'drop'
    ) %>%
    left_join(individual_scan_counts, by = c("individual_a" = "individual", "group", "season_id")) %>%
    rename(scans_individual_a = scans_sampled) %>%
    left_join(individual_scan_counts, by = c("individual_b" = "individual", "group", "season_id")) %>%
    rename(scans_individual_b = scans_sampled) %>%
    mutate(
      total_dyad_scans = scans_individual_a + scans_individual_b - co_occurrence_scans,
      association_rate = co_occurrence_scans / total_dyad_scans
    ) %>%
    arrange(group, season_id, individual_a, individual_b)
  
  return(edgelist)
}
