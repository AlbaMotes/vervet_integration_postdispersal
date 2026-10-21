
#elo rating calculations. It creates both male-male based elos and female-male based elos. 
#It also calculates both at the 1 year post-immigration birthday (target_type="birthday") 
#and at the end of tenure target_type="last_seen"
calculate_elo_ratings <- function(individuals,
                                  conflict_data,
                                  group_presence_matrices,
                                  target_type = c("birthday", "last_seen")) {
  
  library(EloRating)
  target_type <- match.arg(target_type)
  suffix      <- ifelse(target_type == "birthday", "one_year", "end_tenure")
  
  # ── Ensure date columns are Date ──────────────────────────────────────────
  individuals <- individuals %>%
    mutate(DateImmigration1 = as.Date(DateImmigration1),
           LastSeen1        = as.Date(LastSeen1))
  
  conflict_data <- conflict_data %>%
    mutate(date = as.Date(date)) %>%
    arrange(date)
  
  # ── Initialise output columns ─────────────────────────────────────────────
  for (type in c("all", "male_male")) {
    individuals[[paste0("elo_",              suffix, "_", type)]] <- NA_real_
    individuals[[paste0("elo_rating_date_", suffix, "_", type)]] <- as.Date(NA)
    individuals[[paste0("n_conflicts_",     suffix, "_", type)]] <- 0L
  }
  
  # ── Helper: find closest available Elo date to a target ──────────────────
  find_closest_date <- function(target, elo_dates) {
    if (target %in% elo_dates) return(target)
    elo_dates[which.min(abs(elo_dates - target))]
  }
  
  # ── Helper: extract Elo for one individual from one elo.seq object ────────
  extract_elo_for_individual <- function(elo_obj, individual_id,
                                         target_date, individual_conflict_dates) {
    
    # Use truedates slot instead of rownames
    elo_dates <- as.Date(elo_obj$truedates)
    
    if (!individual_id %in% colnames(elo_obj$pmat)) {
      message("Individual ", individual_id, " not found in Elo ratings")
      return(list(rating = NA_real_, date = as.Date(NA)))
    }
    
    # Use last conflict date if target is beyond it, else use target
    if (length(individual_conflict_dates) > 0) {
      last_conflict <- max(individual_conflict_dates, na.rm = TRUE)
      ref_date      <- if (target_date <= last_conflict) target_date else last_conflict
    } else {
      ref_date <- target_date
    }
    
    closest <- find_closest_date(ref_date, elo_dates)
    
    rating <- tryCatch(
      extract_elo(elo_obj, 
                  extractdate = format(closest, "%Y-%m-%d"),
                  IDs         = individual_id, 
                  standardize = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(rating) || length(rating) == 0) return(list(rating = NA_real_, date = as.Date(NA)))
    list(rating = rating, date = closest)
  }
  
  # ── Main loop over individuals ────────────────────────────────────────────
  for (i in seq_len(nrow(individuals))) {
    
    individual_id <- individuals$Code[i]
    immig_group   <- individuals$ImmigrationGp1[i]
    
    if (is.na(immig_group)) next
    
    presence <- group_presence_matrices[[immig_group]]
    if (is.null(presence)) {
      message("No presence data for group: ", immig_group)
      next
    }
    
    presence <- presence %>% mutate(Date = as.Date(Date))
    presence_individuals <- setdiff(colnames(presence), "Date")
    
    # Define time window — prefer DateImmigration1, fall back to
    # DepartureNatalGp when DateImmigration1 is missing
    immig_date <- if (!is.na(individuals$DateImmigration1[i])) {
      individuals$DateImmigration1[i]
    } else if (!is.na(individuals$DepartureNatalGp[i])) {
      as.Date(individuals$DepartureNatalGp[i])
    } else {
      as.Date(NA)
    }

    if (is.na(immig_date)) next

    window_start <- immig_date - 365
    window_end   <- if (target_type == "birthday") immig_date + 365 else individuals$LastSeen1[i]
    target_date  <- window_end

    if (is.na(target_date)) next
    
    # Filter conflicts to this group and time window
    base_conflicts <- conflict_data %>%
      filter(group == immig_group,
             date  >= window_start,
             date  <= window_end,
             winner %in% presence_individuals,
             loser  %in% presence_individuals,
             !is.na(winner_age_cat), winner_age_cat %in% c("AM", "AF"),
             !is.na(loser_age_cat),  loser_age_cat  %in% c("AM", "AF"))
    
    if (nrow(base_conflicts) == 0) {
      message("No conflict data for group: ", immig_group, " for individual: ", individual_id)
      next
    }
    
    # Male-male subset
    mm_conflicts <- base_conflicts %>%
      filter(winner_age_cat == "AM", loser_age_cat == "AM")
    
    # Store conflict counts
    individuals[[paste0("n_conflicts_", suffix, "_all")]][i] <-
      sum(base_conflicts$winner == individual_id | base_conflicts$loser == individual_id)
    individuals[[paste0("n_conflicts_", suffix, "_male_male")]][i] <-
      sum(mm_conflicts$winner == individual_id | mm_conflicts$loser == individual_id)
    
    # ── Run Elo and extract rating for each conflict type ──────────────────
    conflict_sets <- list(all = base_conflicts, male_male = mm_conflicts)
    
    for (type in names(conflict_sets)) {
      
      conflicts <- conflict_sets[[type]]
      
      if (nrow(conflicts) == 0) {
        message("No ", type, " conflicts for group: ", immig_group,
                " for individual: ", individual_id)
        next
      }
      
      # Need at least 2 interactions for elo.seq
      if (nrow(conflicts) < 2) {
        message("Less than 2 ", type, " conflicts for group: ", immig_group,
                " for individual: ", individual_id, " - skipping")
        next
      }
      
      
# Subset presence matrix to only individuals in this conflict set
      conflict_individuals <- unique(c(conflicts$winner, conflicts$loser))
      presence_subset <- presence[, c("Date", intersect(conflict_individuals,
                                                        colnames(presence)))]
      presence_subset <- presence_subset[
        presence_subset$Date >= min(conflicts$date) &
          presence_subset$Date <= max(conflicts$date), 
      ]
      
      elo_obj <- tryCatch(
        elo.seq(winner   = conflicts$winner,
                loser    = conflicts$loser,
                Date     = conflicts$date,
                draw     = NULL,
                presence = presence_subset,
                runcheck = FALSE),
        error = function(e) {
          message("Error calculating Elo (", type, ") for group ", immig_group, ": ", e$message)
          NULL
        }
      )
      
      if (is.null(elo_obj)) next
      
      individual_conflict_dates <- conflicts$date[
        conflicts$winner == individual_id | conflicts$loser == individual_id
      ]
      
      result <- extract_elo_for_individual(elo_obj, individual_id,
                                           target_date, individual_conflict_dates)
      
      individuals[[paste0("elo_",              suffix, "_", type)]][i] <- result$rating
      individuals[[paste0("elo_rating_date_", suffix, "_", type)]][i] <- result$date
    }
  }
  
  return(individuals)
}
