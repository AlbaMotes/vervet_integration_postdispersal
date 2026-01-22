
#elo rating calculations. It creates both male-male based elos and female-male based elos. 
#It also calculates both at the 1 year post-immigration birthday (target_type="birthday") 
#and at the end of tenure target_type="last_seen"

calculate_elo_ratings <- function(individuals, 
                                  conflict_data, 
                                  group_presence_matrices,
                                  target_type = c("birthday", "last_seen")) {
  
  # Validate target_type parameter
  target_type <- match.arg(target_type)
  
  # Convert dates to Date objects if they aren't already
  individuals$DateImmigration1 <- as.Date(individuals$DateImmigration1)
  if(target_type == "last_seen") {
    individuals$LastSeen1 <- as.Date(individuals$LastSeen1)
  }
  conflict_data$date <- as.Date(format(as.Date(conflict_data$date), "%Y-%m-%d"))
  
  # Sort conflict data by date
  conflict_data <- conflict_data[order(conflict_data$date), ]
  
  # Set column name suffix based on target type
  suffix <- ifelse(target_type == "birthday", "one_year", "end_tenure")
  
  # Initialize columns for Elo ratings (all conflicts and male-male only)
  individuals[[paste0("elo_", suffix, "_all")]] <- NA
  individuals[[paste0("elo_", suffix, "_male_male")]] <- NA
  individuals[[paste0("elo_rating_date_", suffix, "_all")]] <- as.Date(NA)
  individuals[[paste0("elo_rating_date_", suffix, "_male_male")]] <- as.Date(NA)
  individuals[[paste0("n_conflicts_", suffix, "_all")]] <- 0
  individuals[[paste0("n_conflicts_", suffix, "_male_male")]] <- 0
  
  # Process each individual
  for(i in 1:nrow(individuals)) {
    # Get individual's ID and immigration group
    individual_id <- individuals$Code[i]
    immig_group <- individuals$ImmigrationGp1[i]
    
    # Get presence data for this group
    presence_data <- group_presence_matrices[[immig_group]]
    
    if(is.null(presence_data)) {
      message(paste("No presence data for group:", immig_group))
      next
    }
    
    # Convert Date column to Date format
    presence_data$Date <- as.Date(format(as.Date(presence_data$Date), "%Y-%m-%d"))
    
    # Get all individuals in the presence data (excluding the Date column)
    presence_individuals <- colnames(presence_data)[-1]
    
    # Make sure presence_data is numeric (except for the Date column)
    for(col in 2:ncol(presence_data)) {
      if(!is.numeric(presence_data[, col])) {
        message(paste("Converting column", col, "to numeric in presence data for group:", immig_group))
        presence_data[, col] <- as.numeric(as.character(presence_data[, col]))
      }
    }
    
    # Get immigration date and define time window
    immig_date <- individuals$DateImmigration1[i]
    one_year_before <- immig_date - 365
    
    # Determine target date and end of window based on target_type
    if(target_type == "birthday") {
      one_year_after <- immig_date + 365
      target_date <- one_year_after
      window_end <- one_year_after
    } else {  # last_seen
      target_date <- individuals$LastSeen1[i]
      window_end <- target_date
    }
    
    # Get conflict data for this group within the time window
    group_conflicts <- conflict_data[
      conflict_data$group == immig_group &
        conflict_data$date >= one_year_before &
        conflict_data$date <= window_end &
        !is.na(conflict_data$winner_age_cat) & conflict_data$winner_age_cat %in% c("AM", "AF") &
        !is.na(conflict_data$loser_age_cat) & conflict_data$loser_age_cat %in% c("AM", "AF"),
    ]
    
    # Skip if no conflict data
    if(nrow(group_conflicts) == 0) {
      message(paste("No conflict data for group:", immig_group, "in the time window for individual:", individual_id))
      next
    }
    
    # Filter out conflicts where either winner or loser is not in presence data
    group_conflicts <- group_conflicts[
      group_conflicts$winner %in% presence_individuals &
        group_conflicts$loser %in% presence_individuals,
    ]
    
    if(nrow(group_conflicts) == 0) {
      message(paste("No valid conflicts after filtering for individual:", individual_id))
      next
    }
    
    # Store number of all conflicts involving this individual
    individuals[[paste0("n_conflicts_", suffix, "_all")]][i] <- sum(
      group_conflicts$winner == individual_id | group_conflicts$loser == individual_id
    )
    
    # Create male-male only conflicts dataset
    group_conflicts_male_male <- group_conflicts[
      !is.na(group_conflicts$winner_age_cat) & group_conflicts$winner_age_cat == "AM" &
        !is.na(group_conflicts$loser_age_cat) & group_conflicts$loser_age_cat == "AM",
    ]
    
    # Store number of male-male conflicts involving this individual
    individuals[[paste0("n_conflicts_", suffix, "_male_male")]][i] <- sum(
      group_conflicts_male_male$winner == individual_id | 
        group_conflicts_male_male$loser == individual_id
    )
    
    # ===== CALCULATE ELO USING ALL CONFLICTS =====
    tryCatch({
      elo_ratings_all <- elo.seq(winner = group_conflicts$winner,
                                 loser = group_conflicts$loser,
                                 Date = group_conflicts$date,
                                 draw = NULL,
                                 presence = presence_data,
                                 runcheck = FALSE)
      
      # Check if individual exists in the Elo ratings
      if(!individual_id %in% colnames(elo_ratings_all$pmat)) {
        message(paste("Individual", individual_id, "not found in Elo ratings (all conflicts) for group", immig_group))
      } else {
        # Skip if invalid date
        if(!is.na(target_date)) {
          # Extract the Elo rating
          tryCatch({
            # Get all dates in the Elo ratings
            elo_dates <- as.Date(rownames(elo_ratings_all$pmat))
            
            # Find the last date this individual appears in the conflict data for this group
            individual_conflict_dates <- group_conflicts$date[
              group_conflicts$winner == individual_id |
                group_conflicts$loser == individual_id
            ]
            
            # Determine which date to use
            if(length(individual_conflict_dates) > 0) {
              last_conflict_date <- max(individual_conflict_dates, na.rm = TRUE)
              
              # Use the target date if it's before or equal to the last conflict date
              # Otherwise use the last conflict date
              if(target_date <= last_conflict_date) {
                # Find the closest date to the target date
                if(target_date %in% elo_dates) {
                  closest_date <- target_date
                } else {
                  closest_idx <- which.min(abs(elo_dates - target_date))
                  closest_date <- elo_dates[closest_idx]
                }
              } else {
                # Use the last conflict date (or closest available in Elo data)
                if(last_conflict_date %in% elo_dates) {
                  closest_date <- last_conflict_date
                } else {
                  closest_idx <- which.min(abs(elo_dates - last_conflict_date))
                  closest_date <- elo_dates[closest_idx]
                }
              }
            } else {
              # If individual never appears in conflict data for this group, use target date
              message(paste("Individual", individual_id, "never appears in conflict data (all) for group", immig_group, "- using target date"))
              if(target_date %in% elo_dates) {
                closest_date <- target_date
              } else {
                closest_idx <- which.min(abs(elo_dates - target_date))
                closest_date <- elo_dates[closest_idx]
              }
            }
            
            # Extract Elo rating for the individual on the closest date
            individual_rating <- extract_elo(elo_ratings_all,
                                             extractdate = closest_date,
                                             IDs = individual_id,
                                             standardize = TRUE)
            
            # Store the rating if it exists
            if(!is.null(individual_rating) && length(individual_rating) > 0) {
              individuals[[paste0("elo_", suffix, "_all")]][i] <- individual_rating
              individuals[[paste0("elo_rating_date_", suffix, "_all")]][i] <- closest_date
            } else {
              message(paste("No rating found for individual", individual_id, "in group", immig_group, "(all conflicts) on closest date"))
            }
          }, error = function(e) {
            message(paste("Error extracting Elo (all conflicts) for individual", individual_id, "in group", immig_group, ":", e$message))
          })
        }
      }
    }, error = function(e) {
      message(paste("Error calculating Elo ratings (all conflicts) for group", immig_group, ":", e$message))
    })
    
    # ===== CALCULATE ELO USING MALE-MALE CONFLICTS ONLY =====
    if(nrow(group_conflicts_male_male) > 0) {
      tryCatch({
        elo_ratings_male_male <- elo.seq(winner = group_conflicts_male_male$winner,
                                         loser = group_conflicts_male_male$loser,
                                         Date = group_conflicts_male_male$date,
                                         draw = NULL,
                                         presence = presence_data,
                                         runcheck = FALSE)
        
        # Check if individual exists in the Elo ratings
        if(!individual_id %in% colnames(elo_ratings_male_male$pmat)) {
          message(paste("Individual", individual_id, "not found in Elo ratings (male-male) for group", immig_group))
        } else {
          # Skip if invalid date
          if(!is.na(target_date)) {
            # Extract the Elo rating
            tryCatch({
              # Get all dates in the Elo ratings
              elo_dates <- as.Date(rownames(elo_ratings_male_male$pmat))
              
              # Find the last date this individual appears in the male-male conflict data
              individual_conflict_dates <- group_conflicts_male_male$date[
                group_conflicts_male_male$winner == individual_id |
                  group_conflicts_male_male$loser == individual_id
              ]
              
              # Determine which date to use
              if(length(individual_conflict_dates) > 0) {
                last_conflict_date <- max(individual_conflict_dates, na.rm = TRUE)
                
                if(target_date <= last_conflict_date) {
                  if(target_date %in% elo_dates) {
                    closest_date <- target_date
                  } else {
                    closest_idx <- which.min(abs(elo_dates - target_date))
                    closest_date <- elo_dates[closest_idx]
                  }
                } else {
                  if(last_conflict_date %in% elo_dates) {
                    closest_date <- last_conflict_date
                  } else {
                    closest_idx <- which.min(abs(elo_dates - last_conflict_date))
                    closest_date <- elo_dates[closest_idx]
                  }
                }
              } else {
                # If individual never appears in male-male conflict data, use target date
                message(paste("Individual", individual_id, "never appears in male-male conflict data for group", immig_group, "- using target date"))
                if(target_date %in% elo_dates) {
                  closest_date <- target_date
                } else {
                  closest_idx <- which.min(abs(elo_dates - target_date))
                  closest_date <- elo_dates[closest_idx]
                }
              }
              
              # Extract Elo rating for the individual on the closest date
              individual_rating <- extract_elo(elo_ratings_male_male,
                                               extractdate = closest_date,
                                               IDs = individual_id,
                                               standardize = TRUE)
              
              # Store the rating if it exists
              if(!is.null(individual_rating) && length(individual_rating) > 0) {
                individuals[[paste0("elo_", suffix, "_male_male")]][i] <- individual_rating
                individuals[[paste0("elo_rating_date_", suffix, "_male_male")]][i] <- closest_date
              } else {
                message(paste("No rating found for individual", individual_id, "in group", immig_group, "(male-male) on closest date"))
              }
            }, error = function(e) {
              message(paste("Error extracting Elo (male-male) for individual", individual_id, "in group", immig_group, ":", e$message))
            })
          }
        }
      }, error = function(e) {
        message(paste("Error calculating Elo ratings (male-male) for group", immig_group, ":", e$message))
      })
    } else {
      message(paste("No male-male conflicts for group:", immig_group, "for individual:", individual_id))
    }
  }
  
  return(individuals)
}