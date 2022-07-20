# Install the necessary packages ----
NPacks <- c("tidyverse", "here", "GCalignR")

# Load required packages ----
pacman::p_load(char = NPacks)

# FUNCTIONS ----
## Function for "manually" moving one peak of one sample one row bellow/above
## on the aligned data frame
move_one_peak <- function(df, Sample, Peak, Dir){
  # Require necessary packages
  require(dplyr)
  require(tidyr)
  
  # Ensure sample is type character
  Sample = Sample %>% as.character()
  
  print(paste("Target peak:", Peak, sep = " "))
  # Move the selected peak for the selected sample
  ## If the movemvent should be toward one row before
  if(Dir == "up"){
    # Define the column where the compound will be moved, based on Dir value
    new_col <- ncol(df %>% select(1:all_of(Peak))) - 1 %>% as.integer()
    
    # Display section of the table for the sample before modification
    print("sample table around target peak before corrective displacement")
    df %>% select(all_of(new_col):(all_of(new_col)+2)) %>%
      filter(rownames(df) == Sample) %>% print()
    
    df[Sample, new_col] <- df[Sample, Peak]
    df[Sample, Peak] <- 0
    
    # Display section of the table for the sample after modification
    print("Sample table around target peak after corrective displacement")
    df %>% select(all_of(new_col):(all_of(new_col)+2)) %>%
      filter(rownames(df) == Sample) %>% print()
  }
  ## If the movemvent should be toward one row after
  if(Dir == "down"){
    new_col <- ncol(df %>% select(1:all_of(Peak))) + 1 %>% as.integer()
    
    # Display section of the table for the sample before modification
    print("sample table around target peak before corrective displacement")
    df %>% select((new_col-2):new_col) %>%
      filter(rownames(df) == Sample) %>% print()
    
    df[Sample, new_col] <- df[Sample, Peak]
    df[Sample, Peak] <- 0
    
    # Display section of the table for the sample after modification
    print("Sample table around target peak after corrective displacement")
    df %>% select((all_of(new_col)-2):all_of(new_col)) %>%
      filter(rownames(df) == Sample) %>% print()
  }
  
  df
}

## Function for moving several peaks of a sample one row bellow/above on the
## aligned data frame. It uses the previously defined function.
move_peaks <- function(df, Sample, peaks_list, movement_dirs){
  df = df
  Sample = Sample
  p_count <- 1
  peaks_movement <- data.frame(Dir = movement_dirs, Peaks = peaks_list)
  for (p in peaks_movement$Peaks) {
    cat('\n')
    print(paste("Peak No.", p_count, sep = " "))
    dir <- peaks_movement %>% filter(Peaks == p) %>% pull(Dir)
    df <- move_one_peak(df, Peak = p, Dir = dir, Sample = Sample)
    p_count <- p_count + 1
    cat('\n')
  }
  print(paste("Finished! The alignment of"
              , p_count-1
              , "peaks was corrected"
              , sep = " "))
  df
}

# Correct miss-alignments ----
# THIS IS PROBABLY NOT NECESSARY FOR THE STANDARDS!!!
## Extract the mean RTs
samples_mean_RT <- samples_RT %>% select(mean_RT)

## Drop mean_RT from the data frames, it will be recalculated on a later step.
### samples
samples_area <- samples_area %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(samples_area)

samples_RT <- samples_RT %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(samples_RT)

### OW
OW_area <- OW_area %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(OW_area)

OW_RT <- OW_RT %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(OW_RT)

## Peaks movement instructions ####
### Create data frames to guide the modifications of each sample
#### samples ####
movements_samples350 <- data.frame(peaks_list = c(paste0("P"
                                                         , c(106, 107, 124)))
                                   , movement_dirs = c('up', 'up', 'up'))
str(movements_samples350)

movements_samples351 <- data.frame(peaks_list = c(paste0("P"
                                                         , c(106, 107, 124)))
                                   , movement_dirs = c('up','up','up'))
str(movements_samples351)

movements_samples352 <- data.frame(peaks_list = c(paste0("P"
                                                         , c(106, 107, 124, 148)))
                                   , movement_dirs = c('up','up','up','up'))
str(movements_samples352)

movements_samples354 <- data.frame(peaks_list =  'P144'
                                   , movement_dirs = 'up')
str(movements_samples354)

#### OW ####
movements_OW328 <- data.frame(peaks_list = c(paste0("P"
                                                    , c(26, 35, 85, 124, 128)))
                              , movement_dirs = c('up', 'up', 'up', 'up', 'up'))
str(movements_OW328)

movements_OW331 <- data.frame(peaks_list = c(paste0("P"
                                                    , c(85, 128)))
                              , movement_dirs = c('up', 'up'))
str(movements_OW331)

movements_OW333 <- data.frame(peaks_list = c(paste0("P"
                                                    , c(52, 128)))
                              , movement_dirs = c('up', 'up'))
str(movements_OW333)

movements_OW337 <- data.frame(peaks_list = c(paste0("P"
                                                    , c(26, 52, 128)))
                              , movement_dirs = c('up', 'up', 'up'))
str(movements_OW337)

movements_OW338 <- data.frame(peaks_list = c(paste0("P"
                                                    , c(52)))
                              , movement_dirs = c('up'))
str(movements_OW338)

movements_OW341 <- data.frame(peaks_list = c(paste0("P"
                                                    , c(52, 124, 128)))
                              , movement_dirs = c('up','up','up'))
str(movements_OW341)

movements_OW345 <- data.frame(peaks_list = c(paste0("P"
                                                    , c(106, 107)))
                              , movement_dirs = c('up', 'up'))
str(movements_OW345)

## Move peaks within samples to correct alignments ####
### The function move_peaks() displaces the specified peaks of a sample 
### by one position up (previous table peak)/down (later table peak)
### It has to be applied to both RT and area data frames
### samples ####
# Sample 350
samples_RT <- move_peaks(samples_RT
                         , Sample = '350'
                         , peaks_list = movements_samples350 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_samples350 %>%
                           pull(movement_dirs))
samples_area <- move_peaks(samples_area
                           , Sample = '350'
                           , peaks_list = movements_samples350 %>%
                             pull(peaks_list)
                           , movement_dirs = movements_samples350 %>%
                             pull(movement_dirs))

# Sample 351
samples_RT <- move_peaks(samples_RT
                         , Sample = '351'
                         , peaks_list = movements_samples351 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_samples351 %>%
                           pull(movement_dirs))
samples_area <- move_peaks(samples_area
                           , Sample = '351'
                           , peaks_list = movements_samples351 %>%
                             pull(peaks_list)
                           , movement_dirs = movements_samples351 %>%
                             pull(movement_dirs))

# Sample 352
samples_RT <- move_peaks(samples_RT
                         , Sample = '352'
                         , peaks_list = movements_samples352 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_samples352 %>%
                           pull(movement_dirs))
samples_area <- move_peaks(samples_area
                           , Sample = '352'
                           , peaks_list = movements_samples352 %>%
                             pull(peaks_list)
                           , movement_dirs = movements_samples352 %>%
                             pull(movement_dirs))

# Sample 354
samples_RT <- move_peaks(samples_RT
                         , Sample = '354'
                         , peaks_list = movements_samples354 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_samples354 %>%
                           pull(movement_dirs))
samples_area <- move_peaks(samples_area
                           , Sample = '354'
                           , peaks_list = movements_samples354 %>%
                             pull(peaks_list)
                           , movement_dirs = movements_samples354 %>%
                             pull(movement_dirs))

### OW ####
# Sample 328
OW_RT <- move_peaks(OW_RT
                    , Sample = '328'
                    , peaks_list = movements_OW328 %>% 
                      pull(peaks_list)
                    , movement_dirs = movements_OW328 %>%
                      pull(movement_dirs))
OW_area <- move_peaks(OW_area
                      , Sample = '328'
                      , peaks_list = movements_OW328 %>% 
                        pull(peaks_list)
                      , movement_dirs = movements_OW328 %>%
                        pull(movement_dirs))

# Sample 331
OW_RT <- move_peaks(OW_RT
                    , Sample = '331'
                    , peaks_list = movements_OW331 %>% 
                      pull(peaks_list)
                    , movement_dirs = movements_OW331 %>%
                      pull(movement_dirs))
OW_area <- move_peaks(OW_area
                      , Sample = '331'
                      , peaks_list = movements_OW331 %>% 
                        pull(peaks_list)
                      , movement_dirs = movements_OW331 %>%
                        pull(movement_dirs))

# Sample 333
OW_RT <- move_peaks(OW_RT
                    , Sample = '333'
                    , peaks_list = movements_OW333 %>% 
                      pull(peaks_list)
                    , movement_dirs = movements_OW333 %>%
                      pull(movement_dirs))
OW_area <- move_peaks(OW_area
                      , Sample = '333'
                      , peaks_list = movements_OW333 %>% 
                        pull(peaks_list)
                      , movement_dirs = movements_OW333 %>%
                        pull(movement_dirs))

# Sample 337
OW_RT <- move_peaks(OW_RT
                    , Sample = '337'
                    , peaks_list = movements_OW337 %>% 
                      pull(peaks_list)
                    , movement_dirs = movements_OW337 %>%
                      pull(movement_dirs))
OW_area <- move_peaks(OW_area
                      , Sample = '337'
                      , peaks_list = movements_OW337 %>% 
                        pull(peaks_list)
                      , movement_dirs = movements_OW337 %>%
                        pull(movement_dirs))

# Sample 338
OW_RT <- move_peaks(OW_RT
                    , Sample = '338'
                    , peaks_list = movements_OW338 %>% 
                      pull(peaks_list)
                    , movement_dirs = movements_OW338 %>%
                      pull(movement_dirs))
OW_area <- move_peaks(OW_area
                      , Sample = '338'
                      , peaks_list = movements_OW338 %>% 
                        pull(peaks_list)
                      , movement_dirs = movements_OW338 %>%
                        pull(movement_dirs))

# Sample 341
OW_RT <- move_peaks(OW_RT
                    , Sample = '341'
                    , peaks_list = movements_OW341 %>%
                      pull(peaks_list)
                    , movement_dirs = movements_OW341 %>%
                      pull(movement_dirs))
OW_area <- move_peaks(OW_area
                      , Sample = '341'
                      , peaks_list = movements_OW341 %>%
                        pull(peaks_list)
                      , movement_dirs = movements_OW341 %>%
                        pull(movement_dirs))


# Sample 345
OW_RT <- move_peaks(OW_RT
                    , Sample = '345'
                    , peaks_list = movements_OW345 %>%
                      pull(peaks_list)
                    , movement_dirs = movements_OW345 %>%
                      pull(movement_dirs))
OW_area <- move_peaks(OW_area
                      , Sample = '345'
                      , peaks_list = movements_OW345 %>%
                        pull(peaks_list)
                      , movement_dirs = movements_OW345 %>%
                        pull(movement_dirs))

print("Alignment was corrected. Check heat map to verify that it is correct.")

# Check the corrected alignment ----
# This is done with the clustered heat map
## Data set with relative abundance of each peak (%)\
### samples
samples_per <- samples_area / rowSums(samples_area) * 100

### OW
OW_per <- OW_area / rowSums(OW_area) * 100

# set.seed(12345)
pdf(here("output", 'corrected-alginment_heatmap.pdf')
    , width = 20, height = 10)
{gplots::heatmap.2(as.matrix(samples_per %>% log1p()), 
                   main = "Clustering of in-hive workers CHC",
                   srtCol = 90,
                   dendrogram = "row",
                   # Rowv = dend,
                   Colv = "NA", # this to make sure the columns are not ordered
                   trace = "none",          
                   # margins =c(5,0.1),      
                   key.xlab = "Relative abundance (%) on Cuticle",
                   #denscol = "grey",
                   #density.info = "density",
                   # RowSideColors = microlab, # to add nice colored strips        
                   col = heatmap_colors)}
{gplots::heatmap.2(as.matrix(OW_per %>% log1p()), 
                   main = "Clustering of out-hive workers CHC",
                   srtCol = 90,
                   dendrogram = "row",
                   # Rowv = dend,
                   Colv = "NA", # this to make sure the columns are not ordered
                   trace = "none",          
                   # margins =c(5,0.1),      
                   key.xlab = "Relative abundance (%) on Cuticle",
                   #denscol = "grey",
                   #density.info = "density",
                   # RowSideColors = microlab, # to add nice colored strips        
                   col = heatmap_colors)}
dev.off()
print("Heat map to verify corrected alignment was exported")

# Export the aligned data for CHC identification ----
## Drop empty columns
### samples
samples_RT <- samples_RT %>% 
  select(-all_of(colnames(samples_RT[colSums(samples_RT) == 0])))

samples_area <- samples_area %>% 
  select(-all_of(colnames(samples_area[colSums(samples_area) == 0])))

samples_mean_RT <- samples_mean_RT %>% t %>% as.data.frame %>% 
  select(all_of(colnames(samples_RT))) %>% t %>% as.data.frame

### OW
OW_RT <- OW_RT %>% 
  select(-all_of(colnames(OW_RT[colSums(OW_RT) == 0])))

OW_area <- OW_area %>% 
  select(-all_of(colnames(OW_area[colSums(OW_area) == 0])))

OW_mean_RT <- OW_mean_RT %>% t %>% as.data.frame %>% 
  select(all_of(colnames(OW_RT))) %>% t %>% as.data.frame
OW_mean_RT

## Recalculate the mean RT
### samples
samples_RT[samples_RT == 0] <- NA
samples_mean_RT <- samples_mean_RT %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = samples_RT %>% colMeans(na.rm = T))
### OW
OW_RT[OW_RT == 0] <- NA
OW_mean_RT <- OW_mean_RT %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = OW_RT %>% colMeans(na.rm = T))

## re-shape data frames before exporting them
### samples
samples_RT <- cbind.data.frame(Peak = row.names(samples_mean_RT)
                               , mean_RT = samples_mean_RT$new_mean_RT
                               # Place samples as columns
                               , samples_RT %>%
                                 t %>% 
                                 as.data.frame %>% 
                                 # Order them in descending order
                                 # , regarding their total abundance
                                 select(all_of(samples_area %>% 
                                                 rowSums() %>% 
                                                 sort(decreasing = T) %>%
                                                 names()))) %>%
  as_tibble()
samples_RT

samples_area <- cbind.data.frame(Peak = row.names(samples_mean_RT)
                                 , mean_RT = samples_mean_RT$new_mean_RT
                                 # Place samples as columns
                                 , samples_area %>%
                                   t %>% 
                                   as.data.frame %>% 
                                   # Order them in descending order
                                   # , regarding their total abundance
                                   select(all_of(samples_area %>% 
                                                   rowSums() %>% 
                                                   sort(decreasing = T) %>%
                                                   names()))) %>%
  as_tibble()
samples_area

### OW
OW_RT <- cbind.data.frame(Peak = row.names(OW_mean_RT)
                          , mean_RT = OW_mean_RT$new_mean_RT
                          # Place samples as columns
                          , OW_RT %>%
                            t %>% 
                            as.data.frame %>% 
                            # Order them in descending order
                            # , regarding their total abundance
                            select(all_of(OW_area %>% 
                                            rowSums() %>% 
                                            sort(decreasing = T) %>%
                                            names()))) %>%
  as_tibble()
OW_RT

OW_area <- cbind.data.frame(Peak = row.names(OW_mean_RT)
                            , mean_RT = OW_mean_RT$new_mean_RT
                            # Place samples as columns
                            , OW_area %>%
                              t %>% 
                              as.data.frame %>% 
                              # Order them in descending order
                              # , regarding their total abundance
                              select(all_of(OW_area %>% 
                                              rowSums() %>% 
                                              sort(decreasing = T) %>%
                                              names()))) %>%
  as_tibble()
OW_area

## Save the data frames 
save(list = c("samples_area", "samples_RT", "OW_area", "OW_RT")
     , file = here("AMelliMelli-CHC-data"
                   , "processed", "aligned_gcms-data.Rdata"))
print("The aligned data frames were exported")

openxlsx::write.xlsx(samples_RT
                     , here("AMelliMelli-CHC-data", "tmp", "samples_RT_table.xlsx"))

openxlsx::write.xlsx(OW_RT
                     , here("AMelliMelli-CHC-data", "tmp", "OW_RT_table.xlsx"))


openxlsx::write.xlsx(STD_RT
                     , here("AMelliMelli-CHC-data"
                            , "tmp", "STD_RT_table.xlsx"))

print("The RT tables for CHC identification was exported")

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "SInf_Script02.txt"))
print("The sessionInfo report was exported. The script 02 finished running")

# ## Detach/unload packages
# lapply(NPacks, unloadNamespace)
# 
# ## Clear environment
# rm(list=ls())


