library(here)

analysis = c("Flying")
gcms_batch = c("1708-2408_2022")
sample_list = c("Samples-Flying.csv")
factors <- c("Subspecies", "Flying")

source(here("scripts", "01_importing_data_automised.R"))

# # Threshold for merging rows step.
# # Skimming the chromatograms no compound seemed to have a separation
# # over 0.157 min between samples
row_merging_threshold <- 0.155

source(here("scripts", "02_align_MS_data_automised.R"))