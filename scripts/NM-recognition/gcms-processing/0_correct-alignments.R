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

# Load the data ----
load(here("data",  "raw"
          , "NM-recognition"
          , "tmp"
          , paste0("uncorrected-alignment_"
                   , "049-271"
                   , ".Rdata")))

# Correct miss-alignments ----
# THIS IS PROBABLY NOT NECESSARY FOR THE STANDARDS!!!
## Extract the mean RTs
Br_Ca_mean_RT <- Br_Ca_RT %>% select(mean_RT)

Br_Ib_mean_RT <- Br_Ib_RT %>% select(mean_RT)

Wu_Ca_mean_RT <- Wu_Ca_RT %>% select(mean_RT)

Wu_Ib_0322_mean_RT <- Wu_Ib_0322_RT %>% select(mean_RT)

Wu_Ib_0422_mean_RT <- Wu_Ib_0422_RT %>% select(mean_RT)

STD_0322_mean_RT <- STD_0322_RT %>% select(mean_RT)

STD_0422_mean_RT <- STD_0422_RT %>% select(mean_RT)

## Drop mean_RT from the data frames, it will be recalculated on a later step.
### Br_Ca
Br_Ca_area <- Br_Ca_area %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(Br_Ca_area)

Br_Ca_RT <- Br_Ca_RT %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(Br_Ca_RT)

### Br_Ib
Br_Ib_area <- Br_Ib_area %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(Br_Ib_area)

Br_Ib_RT <- Br_Ib_RT %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(Br_Ib_RT)

### Wu_Ca
Wu_Ca_area <- Wu_Ca_area %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(Wu_Ca_area)

Wu_Ca_RT <- Wu_Ca_RT %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(Wu_Ca_RT)

### Wu_Ib_0322
Wu_Ib_0322_area <- Wu_Ib_0322_area %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(Wu_Ib_0322_area)

Wu_Ib_0322_RT <- Wu_Ib_0322_RT %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(Wu_Ib_0322_RT)

### Wu_Ib_0422
Wu_Ib_0422_area <- Wu_Ib_0422_area %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(Wu_Ib_0422_area)

Wu_Ib_0422_RT <- Wu_Ib_0422_RT %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(Wu_Ib_0422_RT)

### STD_0322
STD_0322_area <- STD_0322_area %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(STD_0322_area)

STD_0322_RT <- STD_0322_RT %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(STD_0322_RT)

### STD_0422
STD_0422_area <- STD_0422_area %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(STD_0422_area)

STD_0422_RT <- STD_0422_RT %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(STD_0422_RT)

## Peaks movement instructions ####
### Create data frames to guide the modifications of each sample
#### Br_Ca ####
movements_Br_Ca194 <- data.frame(peaks_list = c(paste0("P"
                                                         , c(124)))
                                   , movement_dirs = c('up'))
str(movements_Br_Ca194)

movements_Br_Ca197<- data.frame(peaks_list = c(paste0("P"
                                                       , c(124)))
                                 , movement_dirs = c('up'))
str(movements_Br_Ca197)

#### Wu_Ca ####
movements_Wu_Ca63 <- data.frame(peaks_list = c(paste0("P"
                                                       , c(198)))
                                 , movement_dirs = c('up'))
str(movements_Wu_Ca63)

movements_Wu_Ca64<- data.frame(peaks_list = c(paste0("P"
                                                      , c(198)))
                                , movement_dirs = c('up'))
str(movements_Wu_Ca64)

movements_Wu_Ca87<- data.frame(peaks_list = c(paste0("P"
                                                     , c(198)))
                               , movement_dirs = c('up'))
str(movements_Wu_Ca87)

movements_Wu_Ca94<- data.frame(peaks_list = c(paste0("P"
                                                     , c(198)))
                               , movement_dirs = c('up'))
str(movements_Wu_Ca94)

movements_Wu_Ca166<- data.frame(peaks_list = c(paste0("P"
                                                     , c(198)))
                               , movement_dirs = c('up'))
str(movements_Wu_Ca166)

#### STD_0422 ####
movements_STD_0422_L2204_2 <- data.frame(peaks_list = c(paste0("P"
                                                      , c(8, 87, 105)))
                                , movement_dirs = c('up', 'up', 'up'))
str(movements_STD_0422_L2204_2)

movements_STD_0422_H2204_2<- data.frame(peaks_list = c(paste0("P"
                                                     , c(113, 115, 127, 135
                                                     , 144, 165, 184
                                                     , 187, 189, 192)))
                               , movement_dirs = c('up', 'up', 'up', 'up'
                                                 , 'up', 'up', 'up', 'up'
                                                 , 'up', 'up'))
str(movements_STD_0422_H2204_2)







## Move peaks within samples to correct alignments ####
### The function move_peaks() displaces the specified peaks of a sample 
### by one position up (previous table peak)/down (later table peak)
### It has to be applied to both RT and area data frames
### Br_Ca ####
# Sample 194
Br_Ca_RT <- move_peaks(Br_Ca_RT
                         , Sample = '194'
                         , peaks_list = movements_Br_Ca194 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Br_Ca194 %>%
                           pull(movement_dirs))
Br_Ca_area <- move_peaks(Br_Ca_area
                           , Sample = '194'
                           , peaks_list = movements_Br_Ca194 %>%
                             pull(peaks_list)
                           , movement_dirs = movements_Br_Ca194 %>%
                             pull(movement_dirs))

# Sample 197
Br_Ca_RT <- move_peaks(Br_Ca_RT
                       , Sample = '197'
                       , peaks_list = movements_Br_Ca197 %>%
                         pull(peaks_list)
                       , movement_dirs = movements_Br_Ca197 %>%
                         pull(movement_dirs))
Br_Ca_area <- move_peaks(Br_Ca_area
                         , Sample = '197'
                         , peaks_list = movements_Br_Ca197 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Br_Ca197 %>%
                           pull(movement_dirs))

### Wu_Ca ####
# Sample 63
Wu_Ca_RT <- move_peaks(Wu_Ca_RT
                       , Sample = '63'
                       , peaks_list = movements_Wu_Ca63 %>%
                         pull(peaks_list)
                       , movement_dirs = movements_Wu_Ca63 %>%
                         pull(movement_dirs))
Wu_Ca_area <- move_peaks(Wu_Ca_area
                         , Sample = '63'
                         , peaks_list = movements_Wu_Ca63 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Wu_Ca63 %>%
                           pull(movement_dirs))

# Sample 64
Wu_Ca_RT <- move_peaks(Wu_Ca_RT
                       , Sample = '64'
                       , peaks_list = movements_Wu_Ca64 %>%
                         pull(peaks_list)
                       , movement_dirs = movements_Wu_Ca64 %>%
                         pull(movement_dirs))
Wu_Ca_area <- move_peaks(Wu_Ca_area
                         , Sample = '64'
                         , peaks_list = movements_Wu_Ca64 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Wu_Ca64 %>%
                           pull(movement_dirs))

# Sample 87
Wu_Ca_RT <- move_peaks(Wu_Ca_RT
                       , Sample = '87'
                       , peaks_list = movements_Wu_Ca87 %>%
                         pull(peaks_list)
                       , movement_dirs = movements_Wu_Ca87 %>%
                         pull(movement_dirs))
Wu_Ca_area <- move_peaks(Wu_Ca_area
                         , Sample = '87'
                         , peaks_list = movements_Wu_Ca87 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Wu_Ca87 %>%
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


