#Function to extract proximity scan information of 1m and arm length

extract_proximity_pairs <- function(scans_filtered) {

  proximity_data <- list()

  # Unique identifier per ROW (i.e. per scanned individual within a scan_id),
  # so "already have data for this individual" tracking never accidentally
  # drops other individuals from the same scan_id.
  scans_filtered <- scans_filtered %>% mutate(.row_id = dplyr::row_number())

  data_before_2022 <- scans_filtered %>% filter(year < 2022)
  data_from_2022   <- scans_filtered %>% filter(year >= 2022)

  # Before 2022: 1m priority
  if (nrow(data_before_2022) > 0) {
    processed_rows <- c()

    if (any(!is.na(data_before_2022$ind1m) & data_before_2022$ind1m != "")) {
      ind1m_data <- data_before_2022 %>%
        filter(!is.na(ind1m) & ind1m != "") %>%
        separate_rows(ind1m, sep = ";") %>%
        mutate(interaction_partner = trimws(tolower(ind1m))) %>%
        filter(interaction_partner != "") %>%
        dplyr::select(-ind1m)
      proximity_data <- append(proximity_data, list(ind1m_data))
      processed_rows <- unique(ind1m_data$.row_id)
    }

    remaining <- data_before_2022 %>% filter(!.row_id %in% processed_rows)
    if (nrow(remaining) > 0) {
      dist_1m <- list()
      nnadult_1m <- remaining %>%
        filter(!is.na(distancenna) & distancenna == "1m" & !is.na(nnadult) & nnadult != "") %>%
        separate_rows(nnadult, sep = ";") %>%
        mutate(interaction_partner = trimws(tolower(nnadult))) %>%
        filter(interaction_partner != "") %>% dplyr::select(-nnadult)
      if (nrow(nnadult_1m) > 0) dist_1m <- append(dist_1m, list(nnadult_1m))

      nnjuv_1m <- remaining %>%
        filter(!is.na(distancennj) & distancennj == "1m" & !is.na(nnjuvenile) & nnjuvenile != "") %>%
        separate_rows(nnjuvenile, sep = ";") %>%
        mutate(interaction_partner = trimws(tolower(nnjuvenile))) %>%
        filter(interaction_partner != "") %>% dplyr::select(-nnjuvenile)
      if (nrow(nnjuv_1m) > 0) dist_1m <- append(dist_1m, list(nnjuv_1m))

      if (length(dist_1m) > 0) {
        combined <- bind_rows(dist_1m)
        proximity_data <- append(proximity_data, list(combined))
        processed_rows <- c(processed_rows, unique(combined$.row_id))
      }
    }

    final_remaining <- data_before_2022 %>% filter(!.row_id %in% processed_rows)
    if (nrow(final_remaining) > 0) {
      arm_processed <- c()
      if (any(!is.na(final_remaining$indarmlength) & final_remaining$indarmlength != "")) {
        arm_data <- final_remaining %>%
          filter(!is.na(indarmlength) & indarmlength != "") %>%
          separate_rows(indarmlength, sep = ";") %>%
          mutate(interaction_partner = trimws(tolower(indarmlength))) %>%
          filter(interaction_partner != "") %>% dplyr::select(-indarmlength)
        proximity_data <- append(proximity_data, list(arm_data))
        arm_processed  <- unique(arm_data$.row_id)
      }
      arm_remaining <- final_remaining %>% filter(!.row_id %in% arm_processed)
      if (nrow(arm_remaining) > 0) {
        dist_arm <- list()
        nnadult_arm <- arm_remaining %>%
          filter(!is.na(distancenna) & distancenna == "Arm length" & !is.na(nnadult) & nnadult != "") %>%
          separate_rows(nnadult, sep = ";") %>%
          mutate(interaction_partner = trimws(tolower(nnadult))) %>%
          filter(interaction_partner != "") %>% dplyr::select(-nnadult)
        if (nrow(nnadult_arm) > 0) dist_arm <- append(dist_arm, list(nnadult_arm))

        nnjuv_arm <- arm_remaining %>%
          filter(!is.na(distancennj) & distancennj == "Arm length" & !is.na(nnjuvenile) & nnjuvenile != "") %>%
          separate_rows(nnjuvenile, sep = ";") %>%
          mutate(interaction_partner = trimws(tolower(nnjuvenile))) %>%
          filter(interaction_partner != "") %>% dplyr::select(-nnjuvenile)
        if (nrow(nnjuv_arm) > 0) dist_arm <- append(dist_arm, list(nnjuv_arm))

        if (length(dist_arm) > 0) proximity_data <- append(proximity_data, list(bind_rows(dist_arm)))
      }
    }
  }

  # From 2022: Arm length priority
  if (nrow(data_from_2022) > 0) {
    processed_rows_2022 <- c()

    if (any(!is.na(data_from_2022$indarmlength) & data_from_2022$indarmlength != "")) {
      arm_data_2022 <- data_from_2022 %>%
        filter(!is.na(indarmlength) & indarmlength != "") %>%
        separate_rows(indarmlength, sep = ";") %>%
        mutate(interaction_partner = trimws(tolower(indarmlength))) %>%
        filter(interaction_partner != "") %>% dplyr::select(-indarmlength)
      proximity_data      <- append(proximity_data, list(arm_data_2022))
      processed_rows_2022 <- unique(arm_data_2022$.row_id)
    }

    remaining_2022 <- data_from_2022 %>% filter(!.row_id %in% processed_rows_2022)
    if (nrow(remaining_2022) > 0) {
      dist_arm_2022 <- list()
      nnadult_arm_2022 <- remaining_2022 %>%
        filter(!is.na(distancenna) & distancenna == "Arm length" & !is.na(nnadult) & nnadult != "") %>%
        separate_rows(nnadult, sep = ";") %>%
        mutate(interaction_partner = trimws(tolower(nnadult))) %>%
        filter(interaction_partner != "") %>% dplyr::select(-nnadult)
      if (nrow(nnadult_arm_2022) > 0) dist_arm_2022 <- append(dist_arm_2022, list(nnadult_arm_2022))

      nnjuv_arm_2022 <- remaining_2022 %>%
        filter(!is.na(distancennj) & distancennj == "Arm length" & !is.na(nnjuvenile) & nnjuvenile != "") %>%
        separate_rows(nnjuvenile, sep = ";") %>%
        mutate(interaction_partner = trimws(tolower(nnjuvenile))) %>%
        filter(interaction_partner != "") %>% dplyr::select(-nnjuvenile)
      if (nrow(nnjuv_arm_2022) > 0) dist_arm_2022 <- append(dist_arm_2022, list(nnjuv_arm_2022))

      if (length(dist_arm_2022) > 0) {
        combined_arm_2022   <- bind_rows(dist_arm_2022)
        proximity_data      <- append(proximity_data, list(combined_arm_2022))
        processed_rows_2022 <- c(processed_rows_2022, unique(combined_arm_2022$.row_id))
      }
    }

    final_remaining_2022 <- data_from_2022 %>% filter(!.row_id %in% processed_rows_2022)
    if (nrow(final_remaining_2022) > 0) {
      one_m_processed <- c()
      if (any(!is.na(final_remaining_2022$ind1m) & final_remaining_2022$ind1m != "")) {
        ind1m_2022 <- final_remaining_2022 %>%
          filter(!is.na(ind1m) & ind1m != "") %>%
          separate_rows(ind1m, sep = ";") %>%
          mutate(interaction_partner = trimws(tolower(ind1m))) %>%
          filter(interaction_partner != "") %>% dplyr::select(-ind1m)
        proximity_data  <- append(proximity_data, list(ind1m_2022))
        one_m_processed <- unique(ind1m_2022$.row_id)
      }

      one_m_remaining <- final_remaining_2022 %>% filter(!.row_id %in% one_m_processed)
      if (nrow(one_m_remaining) > 0) {
        dist_1m_2022 <- list()
        nnadult_1m_2022 <- one_m_remaining %>%
          filter(!is.na(distancenna) & distancenna == "1m" & !is.na(nnadult) & nnadult != "") %>%
          separate_rows(nnadult, sep = ";") %>%
          mutate(interaction_partner = trimws(tolower(nnadult))) %>%
          filter(interaction_partner != "") %>% dplyr::select(-nnadult)
        if (nrow(nnadult_1m_2022) > 0) dist_1m_2022 <- append(dist_1m_2022, list(nnadult_1m_2022))

        nnjuv_1m_2022 <- one_m_remaining %>%
          filter(!is.na(distancennj) & distancennj == "1m" & !is.na(nnjuvenile) & nnjuvenile != "") %>%
          separate_rows(nnjuvenile, sep = ";") %>%
          mutate(interaction_partner = trimws(tolower(nnjuvenile))) %>%
          filter(interaction_partner != "") %>% dplyr::select(-nnjuvenile)
        if (nrow(nnjuv_1m_2022) > 0) dist_1m_2022 <- append(dist_1m_2022, list(nnjuv_1m_2022))

        if (length(dist_1m_2022) > 0) proximity_data <- append(proximity_data, list(bind_rows(dist_1m_2022)))
      }
    }
  }

  if (length(proximity_data) == 0) return(NULL)
  bind_rows(proximity_data) %>%
    dplyr::select(-.row_id) %>%
    distinct()
}
