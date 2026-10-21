The repository includes four types of files:
- R scripts
- r formulas
- r environment

R scripts are numbered in the order they should be run inside the bayes_sem folder. File paths need to be adapted to the local user.
Large files are stored in the Zenodo repository: 

all_files <- list.files(pattern = 'scan_part.*\\\\.csv')
scan <- do.call(rbind, lapply(all_files, read.csv))

