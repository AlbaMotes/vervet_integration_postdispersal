The repository includes four types of files:
- R scripts
- csv data files
- r formulas
- r environment

R scripts are numbered in the order they should be run. File paths need to be adapted to the local user.
Large csv files are split into parts that need to be merged before being used as follows (exemplified with the scans data set):

all_files <- list.files(pattern = 'scan_part.*\\\\.csv')
scan <- do.call(rbind, lapply(all_files, read.csv))

