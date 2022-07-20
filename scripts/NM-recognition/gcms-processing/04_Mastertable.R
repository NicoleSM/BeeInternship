# Install pacman package if necessary ----
#install.packages("pacman")

# Load required packages ----
NPacks <- c("tidyverse", "here", "readr")
pacman::p_load(char = NPacks)

# FUNCTIONS ----
## A customized merge function, to ease its usage among the script and within
## other functions
my_merge <- function(x, y, merge_by) {
  merge(x = x, y = y
        , by = all_of(merge_by)
        , all = T
        , sort = F)
}

## Function to recalculate the RI in the first master table
mt_new_ri <- function(master.table) {
  # Which compounds are duplicated/misaligned between group tables
  # inside the master table?
  dup_comps <- master.table %>% filter(duplicated(Compound)) %>%
    select(RI, Compound)
  str(dup_comps)
  
  # Prepare to recalculate new RIs
  dup_comps2 <- master.table %>% distinct(Compound, .keep_all = T) %>%
    filter(Compound %in% dup_comps$Compound) %>%
    select(RI, Compound)
  dup_comps <- rbind.data.frame(dup_comps, dup_comps2) %>% arrange(RI)
  dup_comps <- dup_comps %>%
    mutate(new.RIs = rep(NA, nrow(dup_comps))) %>%
    as.data.frame()
  rm(dup_comps2)
  
  # Loop to recalculate RIs of repeated compounds allowing proper merge
  # It iterates on each unique compound name in dup_comps data frame
  Lap = 1
  for (i in unique(dup_comps$Compound)) {
    # Report current lap compound
    print(paste("Lap", Lap, ":", i, sep = " "))
    # Store peaks with the compound in a new temporal data frame
    peaks <- dup_comps %>% filter(Compound == i)
    # Store the RIs in descending order in a new temporal data frame
    desc_RIs <- peaks %>% arrange(desc(RI)) %>% select(RI)
    # The desc_RIs data frame is transformed into a string
    desc_RIs <- desc_RIs["RI"]
    
    # desc_RIs is used to determine the order of iterations
    # for the recalculation of the RIs inside a for loop
    for (n in desc_RIs[, "RI"]) {
      # Report RI determining the current loop iteration
      print(paste("RI:", n, sep = " "))
      # Store peaks that do not differ in their RIs by more than 
      # 2 in a new temporal data frame
      same_peak <- peaks %>%
        filter((n-peaks$RI) <= 2 & (n-peaks$RI) >= 0)
      same_peak
      # Calculate the new RI for the peak, if it has not been done
      # in a previous iteration
      if(is.na(same_peak$new.RIs[][same_peak$RI == n])){
        # New RI is the average of the RIs of the peaks with the same
        # compound name (curren iteration compound)
        # , for which the alignment has to be corrected
        new_RI <- same_peak$RI %>%  mean() %>%
          round(digits = 0) %>%
          rep(nrow(same_peak))
        new_RI
        peaks$new.RIs[][peaks$RI %in% same_peak$RI] <- new_RI
        print(paste("New RI :"
                    , min(new_RI)
                    , sep = " "))
      } else {
        print(paste("New RI is the same as previous ("
                    , min(new_RI), ")"))
      }
    }
    dup_comps$new.RIs[][dup_comps$RI[] %in% peaks$RI] <- peaks$new.RIs
    Lap = Lap + 1
  }
  
  # Replace the RIs with the new RIs
  master.table$RI[][master.table$RI[] %in% dup_comps$RI] <-
    dup_comps$new.RIs[]
  print("All new RIs have been stored in the Master table")
  master.table
}

## Functions to fuse peaks
### Function to fuse the indicated peaks
fuse_peaks <- function(master.table, peaks_to_fuse){
  require(dplyr)
  
  peaks_sum <- master.table %>% 
    filter(Peak %in% peaks_to_fuse) %>% 
    select(!Peak:RI) %>% 
    colSums(na.rm = T) %>% 
    t() %>% 
    as.data.frame()
  
  empty_peaks <- matrix(NA
                        , ncol = length(peaks_sum)
                        , nrow = length(peaks_to_fuse) - 1) %>% 
    as.data.frame()
  colnames(empty_peaks) <- colnames(peaks_sum)
  
  peaks_sum <- rbind(peaks_sum
                     , empty_peaks)
  
  new_RI <- master.table %>% 
    filter(Peak %in% peaks_to_fuse) %>% 
    select(RI) %>% 
    colMeans() %>% 
    t() %>% 
    as.data.frame()
  new_RI <- data.frame(RI = rep(new_RI %>% 
                                   pull(RI)
                                 , length(peaks_to_fuse)))
  
  peaks_sum <- cbind(master.table %>% 
                       filter(Peak %in% peaks_to_fuse) %>% 
                       select(Peak:Mod.position)
                     , new_RI
                     , peaks_sum) %>% 
    as_tibble()
  
  peaks_sum$RI <- peaks_sum$RI %>% as.integer()
  
  master.table <- rows_update(master.table
                         , peaks_sum
                         , by = 'Peak')
  master.table
}

### Function for performing all the fusion of peaks indicated by a list.
### It uses the previously defined function to perform each fusion.
fuse_all_peaks <- function(master.table, fusion.list){
  f_count <- 1
  for(i in fusion.list){
    cat('\n')
    print(paste("Fusion No.", f_count, sep = " "))
    print(paste(length(i)
                ,"peaks to be fused:"
                , paste(i, collapse = " ")
                , sep = " "))
    cat('\n')
    print("Peaks before fusion")
    master.table %>% filter(Peak %in% i) %>% print()
    master.table <- fuse_peaks(master.table = master.table, peaks_to_fuse = i)
    print("Peaks after fusion")
    master.table %>% filter(Peak %in% i) %>% print()
    f_count <- f_count + 1
    cat('\n')
  }
  print(paste("Finished!"
              , f_count-1
              , "fusions were performed"
              , sep = " "))
  master.table
}

# Load group tables ----
load(here("AMelliMelli-CHC-data", "processed", "group_tables.Rdata"))

# First MT ----
comps_vars <- IW_table %>% 
  select(Compound:RI) %>% 
  colnames() # c("Compound", "Class", "Chain.length", "RI")
comps_vars

master.table <- my_merge(IW_table
         , OW_table
         , merge_by = comps_vars) %>% arrange(RI)
str(master.table)

# Export first master table
write_csv(master.table, here("AMelliMelli-CHC-data"
                             , "tmp"
                             , "master_table_first.csv"))

## RI  recalculation ####
master.table <- mt_new_ri(master.table)
str(master.table)

# Export the master table with recalculated RIs
write_csv(master.table, here("AMelliMelli-CHC-data"
                             , "tmp"
                             , "master_table_new-RI.csv"))

# New group tables ----
IW_table <- master.table %>% select(all_of(colnames(IW_table)))
IW_table[is.na(IW_table)] <- 0

OW_table <- master.table %>% select(all_of(colnames(OW_table)))
OW_table[is.na(OW_table)] <- 0

IW_table <- IW_table %>% 
  filter(IW_table %>% 
           select(-all_of(comps_vars)) %>% 
           rowSums() > 0.0)
str(IW_table)

OW_table <- OW_table %>% 
  filter(OW_table %>% 
           select(-all_of(comps_vars)) %>% 
           rowSums() > 0.0)
str(OW_table)

# Second MT ----
# Merge the data frames to build the master table
master.table <- my_merge(IW_table
                         , OW_table
                         , merge_by = comps_vars) %>% 
  arrange(RI) 

master.table <- master.table %>% 
  mutate(Peak = paste0("P", 1:nrow(master.table))) %>% 
  select(Peak, all_of(comps_vars), everything())
str(master.table)

# Fuse peaks ----
## Create a list of fusion operations to be performed
## Each entry of the list corresponds to the names of two or more peaks that
## have to be fused together across the data set.
fusion_list <- list(c(paste0("P"
                             , c(1, 2)))
                    , c(paste0("P"
                               , c(22, 23, 24)))
                    , c(paste0("P"
                               , c(27, 28)))
                    , c(paste0("P"
                               , c(30, 31)))
                    , c(paste0("P"
                               , c(44, 45)))
                    , c(paste0("P"
                               , c(48, 49, 50, 51)))
                    , c(paste0("P"
                               , c(52, 53)))
                    , c(paste0("P"
                               , c(57, 58)))
                    , c(paste0("P"
                               , c(59, 60))))
fusion_list

## Fuse the corresponding peaks using the fuse_all_areas/RTs functions
master.table <- fuse_all_peaks(master.table, fusion_list)
master.table

# Export the master table with fused peaks
write_csv(master.table, here("AMelliMelli-CHC-data"
                             , "tmp"
                             , "master_table_fused-peaks.csv"))

## Drop empty rows
master.table <- master.table %>% 
  filter(master.table %>% 
           select(!Peak:RI) %>% 
           rowSums(na.rm = T) > 0)

# Export data frames ----
master.table <- master.table %>% 
  mutate(Peak = paste0("P", 1:nrow(master.table)))
row.names(master.table) <- master.table$Peak
master.table <- master.table %>% select(-Peak)

master.daten <- master.table %>% 
  select(!all_of(comps_vars)) %>% 
  t() %>% 
  as.data.frame()
master.daten[is.na(master.daten)] <- 0
str(master.daten)

master.Comps <- master.table %>% 
  select(all_of(comps_vars))
master.Comps$RI <- master.Comps$RI %>% as.integer()
str(master.Comps)

save(list = c("master.daten"
              , "master.Comps"
              , "grouping_info")
     , file = here("AMelliMelli-CHC-data"
                   , "processed"
                   , "data-frames.Rdata"))

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "sInf_Script04.txt"))

# ## Detach/unload packages
# lapply(NPacks, unloadNamespace)
# 
# ## Clear environment
# rm(list=ls())
