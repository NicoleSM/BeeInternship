library(here)

analysis = c("queen-less")
sample_list = c("Queenless_hive_table.csv")
factors <- c("Task", "Date")

source(here("scripts", "01_importing_data_automised.R"))
