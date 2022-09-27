# Install the necessary packages ----
NPacks <- c("tidyverse", "here", "GCalignR")

# Load required packages ----
pacman::p_load(char = NPacks)

# FUNCTIONS ----
# mean_RT ####
## A function to obtain mean RT of group dataframes within the master group lis
mean_RT_df <- function(df){
  df_mean_RT <- df %>% select(mean_RT)
  df_mean_RT
  
}
## Function to drop mean_RT from dataframes
drop_mean_RT_df <- function(df){
  df_dropped <- df %>% 
    select(-mean_RT) %>% 
    t %>% 
    as.data.frame()
  str(df_dropped)
  df_dropped
}


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
          , analysis
          , "tmp"
          , paste0("uncorrected-alignment_"
                   , gcms_batch
                   , ".Rdata")))

# Correct miss-alignments ----
# THIS IS PROBABLY NOT NECESSARY FOR THE STANDARDS!!!
## Extract the mean RTs
mg_list_mean_RT <-  lapply(mg_list_RT, mean_RT_df)
print(mg_list_mean_RT)

STD_mean_RT <- STD_RT %>% select(mean_RT)


## Drop mean_RT from the data frames, it will be recalculated on a later step.
### mg_list

mg_list_area <- lapply(mg_list_area, drop_mean_RT_df)

mg_list_RT <- lapply(mg_list_RT, drop_mean_RT_df)

### STD
STD_RT <- STD_RT %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(STD_RT)


## Peaks movement instructions ####
### Create data frames to guide the modifications of each sample
#### Ca_FALSE ####
movements_Ca_FALSE20 <- data.frame(peaks_list = c(paste0("P"
                                                         , c(140, 161)))
                                   , movement_dirs = c('up', 'up'))
str(movements_Ca_FALSE20)

#### Ca_TRUE ####
movements_Ca_TRUE10 <- data.frame(peaks_list = c(paste0("P"
                                                       , c(125)))
                                 , movement_dirs = c('up'))
str(movements_Ca_TRUE10)

movements_Ca_TRUE07<- data.frame(peaks_list = c(paste0("P"
                                                      , c(131, 86, 85, 106)))
                                , movement_dirs = c('up', 'down', 'down', 'up'))
str(movements_Ca_TRUE07)

movements_Ca_TRUE13<- data.frame(peaks_list = c(paste0("P"
                                                       , c(85)))
                                 , movement_dirs = c('down'))
str(movements_Ca_TRUE13)

movements_Ca_TRUE14<- data.frame(peaks_list = c(paste0("P"
                                                       , c(85)))
                                 , movement_dirs = c('down'))
str(movements_Ca_TRUE14)

movements_Ca_TRUE15<- data.frame(peaks_list = c(paste0("P"
                                                       , c(85)))
                                 , movement_dirs = c('down'))
str(movements_Ca_TRUE15)

movements_Ca_TRUE03<- data.frame(peaks_list = c(paste0("P"
                                                       , c(85)))
                                 , movement_dirs = c('down'))
str(movements_Ca_TRUE03)

movements_Ca_TRUE12<- data.frame(peaks_list = c(paste0("P"
                                                       , c(86, 85)))
                                 , movement_dirs = c('down','down'))
str(movements_Ca_TRUE12)

#### Ib_FALSE ####
movements_Ib_FALSE21 <- data.frame(peaks_list = c(paste0("P"
                                                        , c(173, 187, 164)))
                                  , movement_dirs = c('down','up', 'up'))
str(movements_Ib_FALSE21)

movements_Ib_FALSE23<- data.frame(peaks_list = c(paste0("P"
                                                       , c(173)))
                                 , movement_dirs = c('down'))
str(movements_Ib_FALSE23)

movements_Ib_FALSE24<- data.frame(peaks_list = c(paste0("P"
                                                        , c(173, 164)))
                                  , movement_dirs = c('down', 'up'))
str(movements_Ib_FALSE24)

movements_Ib_FALSE25<- data.frame(peaks_list = c(paste0("P"
                                                        , c(173, 187, 164)))
                                  , movement_dirs = c('down','up', 'up'))
str(movements_Ib_FALSE25)

movements_Ib_FALSE26<- data.frame(peaks_list = c(paste0("P"
                                                        , c(173, 187, 164)))
                                  , movement_dirs = c('down','up', 'up'))
str(movements_Ib_FALSE26)

movements_Ib_FALSE31<- data.frame(peaks_list = c(paste0("P"
                                                        , c(173)))
                                  , movement_dirs = c('down'))
str(movements_Ib_FALSE31)

movements_Ib_FALSE36<- data.frame(peaks_list = c(paste0("P"
                                                        , c(173)))
                                  , movement_dirs = c('down'))
str(movements_Ib_FALSE36)

movements_Ib_FALSE32<- data.frame(peaks_list = c(paste0("P"
                                                        , c(256)))
                                  , movement_dirs = c('up'))
str(movements_Ib_FALSE32)

#### Ib_TRUE ####
movements_Ib_TRUE27 <- data.frame(peaks_list = c(paste0("P"
                                                         , c(99, 164)))
                                   , movement_dirs = c('up', 'down'))
str(movements_Ib_TRUE27)

movements_Ib_TRUE28 <- data.frame(peaks_list = c(paste0("P"
                                                        , c(99, 209, 164)))
                                  , movement_dirs = c('up', 'up', 'down'))
str(movements_Ib_TRUE28)

movements_Ib_TRUE29 <- data.frame(peaks_list = c(paste0("P"
                                                        , c(99, 164)))
                                  , movement_dirs = c('up', 'down'))
str(movements_Ib_TRUE29)

movements_Ib_TRUE39 <- data.frame(peaks_list = c(paste0("P"
                                                        , c(99, 138, 164)))
                                  , movement_dirs = c('up', 'up', 'down'))
str(movements_Ib_TRUE39)

movements_Ib_TRUE38 <- data.frame(peaks_list = c(paste0("P"
                                                        , c(139, 138, 144)))
                                  , movement_dirs = c('up', 'up', 'up'))
str(movements_Ib_TRUE38)

movements_Ib_TRUE30 <- data.frame(peaks_list = c(paste0("P"
                                                        , c(163, 138, 164)))
                                  , movement_dirs = c('down', 'up', 'down'))
str(movements_Ib_TRUE30)

movements_Ib_TRUE22 <- data.frame(peaks_list = c(paste0("P"
                                                        , c(144, 126, 150, 206)))
                                  , movement_dirs = c('up', 'up', 'down', 'down'))
str(movements_Ib_TRUE22)

movements_Ib_TRUE40 <- data.frame(peaks_list = c(paste0("P"
                                                        , c(144, 209, 163)))
                                  , movement_dirs = c('up', 'down', 'down'))
str(movements_Ib_TRUE40)

movements_Ib_TRUE33 <- data.frame(peaks_list = c(paste0("P"
                                                        , c(164)))
                                  , movement_dirs = c('down'))
str(movements_Ib_TRUE33)

movements_Ib_TRUE37 <- data.frame(peaks_list = c(paste0("P"
                                                        , c(164)))
                                  , movement_dirs = c('down'))
str(movements_Ib_TRUE37)





## Move peaks within samples to correct alignments ####
### The function move_peaks() displaces the specified peaks of a sample 
### by one position up (previous table peak)/down (later table peak)
### It has to be applied to both RT and area data frames
### Ca_FALSE ####
# Sample 20
mg_list_RT[["Ca_FALSE"]] <- move_peaks(mg_list_RT[["Ca_FALSE"]]
                         , Sample = 'Fly_020'
                         , peaks_list = movements_Ca_FALSE20 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Ca_FALSE20 %>%
                           pull(movement_dirs))
mg_list_area[["Ca_FALSE"]] <- move_peaks(mg_list_area[["Ca_FALSE"]]
                           , Sample = 'Fly_020'
                           , peaks_list = movements_Ca_FALSE20 %>%
                             pull(peaks_list)
                           , movement_dirs = movements_Ca_FALSE20 %>%
                             pull(movement_dirs))
### Ca_TRUE ####
# Sample 07
mg_list_RT[["Ca_TRUE"]] <- move_peaks(mg_list_RT[["Ca_TRUE"]]
                                       , Sample = 'Fly_007'
                                       , peaks_list = movements_Ca_TRUE07 %>%
                                         pull(peaks_list)
                                       , movement_dirs = movements_Ca_TRUE07 %>%
                                         pull(movement_dirs))
mg_list_area[["Ca_TRUE"]] <- move_peaks(mg_list_area[["Ca_TRUE"]]
                                         , Sample = 'Fly_007'
                                         , peaks_list = movements_Ca_TRUE07 %>%
                                           pull(peaks_list)
                                         , movement_dirs = movements_Ca_TRUE07 %>%
                                           pull(movement_dirs))
# Sample 10
mg_list_RT[["Ca_TRUE"]] <- move_peaks(mg_list_RT[["Ca_TRUE"]]
                                      , Sample = 'Fly_010'
                                      , peaks_list = movements_Ca_TRUE10 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ca_TRUE10 %>%
                                        pull(movement_dirs))
mg_list_area[["Ca_TRUE"]] <- move_peaks(mg_list_area[["Ca_TRUE"]]
                                        , Sample = 'Fly_010'
                                        , peaks_list = movements_Ca_TRUE10 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ca_TRUE10 %>%
                                          pull(movement_dirs))
# Sample 13
mg_list_RT[["Ca_TRUE"]] <- move_peaks(mg_list_RT[["Ca_TRUE"]]
                                      , Sample = 'Fly_013'
                                      , peaks_list = movements_Ca_TRUE13 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ca_TRUE13 %>%
                                        pull(movement_dirs))
mg_list_area[["Ca_TRUE"]] <- move_peaks(mg_list_area[["Ca_TRUE"]]
                                        , Sample = 'Fly_013'
                                        , peaks_list = movements_Ca_TRUE13 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ca_TRUE13 %>%
                                          pull(movement_dirs))
# Sample 14
mg_list_RT[["Ca_TRUE"]] <- move_peaks(mg_list_RT[["Ca_TRUE"]]
                                      , Sample = 'Fly_014'
                                      , peaks_list = movements_Ca_TRUE14 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ca_TRUE14 %>%
                                        pull(movement_dirs))
mg_list_area[["Ca_TRUE"]] <- move_peaks(mg_list_area[["Ca_TRUE"]]
                                        , Sample = 'Fly_014'
                                        , peaks_list = movements_Ca_TRUE14 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ca_TRUE14 %>%
                                          pull(movement_dirs))
# Sample 15
mg_list_RT[["Ca_TRUE"]] <- move_peaks(mg_list_RT[["Ca_TRUE"]]
                                      , Sample = 'Fly_015'
                                      , peaks_list = movements_Ca_TRUE15 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ca_TRUE15 %>%
                                        pull(movement_dirs))
mg_list_area[["Ca_TRUE"]] <- move_peaks(mg_list_area[["Ca_TRUE"]]
                                        , Sample = 'Fly_015'
                                        , peaks_list = movements_Ca_TRUE15 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ca_TRUE15 %>%
                                          pull(movement_dirs))
# Sample 03
mg_list_RT[["Ca_TRUE"]] <- move_peaks(mg_list_RT[["Ca_TRUE"]]
                                      , Sample = 'Fly_003'
                                      , peaks_list = movements_Ca_TRUE03 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ca_TRUE03 %>%
                                        pull(movement_dirs))
mg_list_area[["Ca_TRUE"]] <- move_peaks(mg_list_area[["Ca_TRUE"]]
                                        , Sample = 'Fly_003'
                                        , peaks_list = movements_Ca_TRUE03 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ca_TRUE03 %>%
                                          pull(movement_dirs))
# Sample 12
mg_list_RT[["Ca_TRUE"]] <- move_peaks(mg_list_RT[["Ca_TRUE"]]
                                      , Sample = 'Fly_012'
                                      , peaks_list = movements_Ca_TRUE12 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ca_TRUE12 %>%
                                        pull(movement_dirs))
mg_list_area[["Ca_TRUE"]] <- move_peaks(mg_list_area[["Ca_TRUE"]]
                                        , Sample = 'Fly_012'
                                        , peaks_list = movements_Ca_TRUE12 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ca_TRUE12 %>%
                                          pull(movement_dirs))
### Ib_FALSE ####
# Sample 21
mg_list_RT[["Ib_FALSE"]] <- move_peaks(mg_list_RT[["Ib_FALSE"]]
                                      , Sample = 'Fly_021'
                                      , peaks_list = movements_Ib_FALSE21 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ib_FALSE21 %>%
                                        pull(movement_dirs))
mg_list_area[["Ib_FALSE"]] <- move_peaks(mg_list_area[["Ib_FALSE"]]
                                        , Sample = 'Fly_021'
                                        , peaks_list = movements_Ib_FALSE21 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ib_FALSE21 %>%
                                          pull(movement_dirs))

# Sample 23
mg_list_RT[["Ib_FALSE"]] <- move_peaks(mg_list_RT[["Ib_FALSE"]]
                                       , Sample = 'Fly_023'
                                       , peaks_list = movements_Ib_FALSE23 %>%
                                         pull(peaks_list)
                                       , movement_dirs = movements_Ib_FALSE23 %>%
                                         pull(movement_dirs))
mg_list_area[["Ib_FALSE"]] <- move_peaks(mg_list_area[["Ib_FALSE"]]
                                         , Sample = 'Fly_023'
                                         , peaks_list = movements_Ib_FALSE23 %>%
                                           pull(peaks_list)
                                         , movement_dirs = movements_Ib_FALSE23 %>%
                                           pull(movement_dirs))

# Sample 24
mg_list_RT[["Ib_FALSE"]] <- move_peaks(mg_list_RT[["Ib_FALSE"]]
                                       , Sample = 'Fly_024'
                                       , peaks_list = movements_Ib_FALSE24 %>%
                                         pull(peaks_list)
                                       , movement_dirs = movements_Ib_FALSE24 %>%
                                         pull(movement_dirs))
mg_list_area[["Ib_FALSE"]] <- move_peaks(mg_list_area[["Ib_FALSE"]]
                                         , Sample = 'Fly_024'
                                         , peaks_list = movements_Ib_FALSE24 %>%
                                           pull(peaks_list)
                                         , movement_dirs = movements_Ib_FALSE24 %>%
                                           pull(movement_dirs))

# Sample 25
mg_list_RT[["Ib_FALSE"]] <- move_peaks(mg_list_RT[["Ib_FALSE"]]
                                       , Sample = 'Fly_025'
                                       , peaks_list = movements_Ib_FALSE25 %>%
                                         pull(peaks_list)
                                       , movement_dirs = movements_Ib_FALSE25 %>%
                                         pull(movement_dirs))
mg_list_area[["Ib_FALSE"]] <- move_peaks(mg_list_area[["Ib_FALSE"]]
                                         , Sample = 'Fly_025'
                                         , peaks_list = movements_Ib_FALSE25 %>%
                                           pull(peaks_list)
                                         , movement_dirs = movements_Ib_FALSE25 %>%
                                           pull(movement_dirs))
# Sample 26
mg_list_RT[["Ib_FALSE"]] <- move_peaks(mg_list_RT[["Ib_FALSE"]]
                                       , Sample = 'Fly_026'
                                       , peaks_list = movements_Ib_FALSE26 %>%
                                         pull(peaks_list)
                                       , movement_dirs = movements_Ib_FALSE26 %>%
                                         pull(movement_dirs))
mg_list_area[["Ib_FALSE"]] <- move_peaks(mg_list_area[["Ib_FALSE"]]
                                         , Sample = 'Fly_026'
                                         , peaks_list = movements_Ib_FALSE26 %>%
                                           pull(peaks_list)
                                         , movement_dirs = movements_Ib_FALSE26 %>%
                                           pull(movement_dirs))

# Sample 31
mg_list_RT[["Ib_FALSE"]] <- move_peaks(mg_list_RT[["Ib_FALSE"]]
                                       , Sample = 'Fly_031'
                                       , peaks_list = movements_Ib_FALSE31 %>%
                                         pull(peaks_list)
                                       , movement_dirs = movements_Ib_FALSE31 %>%
                                         pull(movement_dirs))
mg_list_area[["Ib_FALSE"]] <- move_peaks(mg_list_area[["Ib_FALSE"]]
                                         , Sample = 'Fly_031'
                                         , peaks_list = movements_Ib_FALSE31 %>%
                                           pull(peaks_list)
                                         , movement_dirs = movements_Ib_FALSE31 %>%
                                           pull(movement_dirs))

# Sample 32
mg_list_RT[["Ib_FALSE"]] <- move_peaks(mg_list_RT[["Ib_FALSE"]]
                                       , Sample = 'Fly_032'
                                       , peaks_list = movements_Ib_FALSE32 %>%
                                         pull(peaks_list)
                                       , movement_dirs = movements_Ib_FALSE32 %>%
                                         pull(movement_dirs))
mg_list_area[["Ib_FALSE"]] <- move_peaks(mg_list_area[["Ib_FALSE"]]
                                         , Sample = 'Fly_032'
                                         , peaks_list = movements_Ib_FALSE32 %>%
                                           pull(peaks_list)
                                         , movement_dirs = movements_Ib_FALSE32 %>%
                                           pull(movement_dirs))

# Sample 36
mg_list_RT[["Ib_FALSE"]] <- move_peaks(mg_list_RT[["Ib_FALSE"]]
                                       , Sample = 'Fly_036'
                                       , peaks_list = movements_Ib_FALSE36 %>%
                                         pull(peaks_list)
                                       , movement_dirs = movements_Ib_FALSE36 %>%
                                         pull(movement_dirs))
mg_list_area[["Ib_FALSE"]] <- move_peaks(mg_list_area[["Ib_FALSE"]]
                                         , Sample = 'Fly_036'
                                         , peaks_list = movements_Ib_FALSE36 %>%
                                           pull(peaks_list)
                                         , movement_dirs = movements_Ib_FALSE36 %>%
                                           pull(movement_dirs))
### Ib_TRUE ####
# Sample 22
mg_list_RT[["Ib_TRUE"]] <- move_peaks(mg_list_RT[["Ib_TRUE"]]
                                      , Sample = 'Fly_022'
                                      , peaks_list = movements_Ib_TRUE22 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ib_TRUE22 %>%
                                        pull(movement_dirs))
mg_list_area[["Ib_TRUE"]] <- move_peaks(mg_list_area[["Ib_TRUE"]]
                                        , Sample = 'Fly_022'
                                        , peaks_list = movements_Ib_TRUE22 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ib_TRUE22 %>%
                                          pull(movement_dirs))
# Sample 27
mg_list_RT[["Ib_TRUE"]] <- move_peaks(mg_list_RT[["Ib_TRUE"]]
                                      , Sample = 'Fly_027'
                                      , peaks_list = movements_Ib_TRUE27 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ib_TRUE27 %>%
                                        pull(movement_dirs))
mg_list_area[["Ib_TRUE"]] <- move_peaks(mg_list_area[["Ib_TRUE"]]
                                        , Sample = 'Fly_027'
                                        , peaks_list = movements_Ib_TRUE27 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ib_TRUE27 %>%
                                          pull(movement_dirs))
# Sample 28
mg_list_RT[["Ib_TRUE"]] <- move_peaks(mg_list_RT[["Ib_TRUE"]]
                                      , Sample = 'Fly_028'
                                      , peaks_list = movements_Ib_TRUE28 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ib_TRUE28 %>%
                                        pull(movement_dirs))
mg_list_area[["Ib_TRUE"]] <- move_peaks(mg_list_area[["Ib_TRUE"]]
                                        , Sample = 'Fly_028'
                                        , peaks_list = movements_Ib_TRUE28 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ib_TRUE28 %>%
                                          pull(movement_dirs))
# Sample 29
mg_list_RT[["Ib_TRUE"]] <- move_peaks(mg_list_RT[["Ib_TRUE"]]
                                      , Sample = 'Fly_029'
                                      , peaks_list = movements_Ib_TRUE29 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ib_TRUE29 %>%
                                        pull(movement_dirs))
mg_list_area[["Ib_TRUE"]] <- move_peaks(mg_list_area[["Ib_TRUE"]]
                                        , Sample = 'Fly_029'
                                        , peaks_list = movements_Ib_TRUE29 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ib_TRUE29 %>%
                                          pull(movement_dirs))
# Sample 30
mg_list_RT[["Ib_TRUE"]] <- move_peaks(mg_list_RT[["Ib_TRUE"]]
                                      , Sample = 'Fly_030'
                                      , peaks_list = movements_Ib_TRUE30 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ib_TRUE30 %>%
                                        pull(movement_dirs))
mg_list_area[["Ib_TRUE"]] <- move_peaks(mg_list_area[["Ib_TRUE"]]
                                        , Sample = 'Fly_030'
                                        , peaks_list = movements_Ib_TRUE30 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ib_TRUE30 %>%
                                          pull(movement_dirs))
# Sample 38
mg_list_RT[["Ib_TRUE"]] <- move_peaks(mg_list_RT[["Ib_TRUE"]]
                                      , Sample = 'Fly_038'
                                      , peaks_list = movements_Ib_TRUE38 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ib_TRUE38 %>%
                                        pull(movement_dirs))
mg_list_area[["Ib_TRUE"]] <- move_peaks(mg_list_area[["Ib_TRUE"]]
                                        , Sample = 'Fly_038'
                                        , peaks_list = movements_Ib_TRUE38 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ib_TRUE38 %>%
                                          pull(movement_dirs))
# Sample 39
mg_list_RT[["Ib_TRUE"]] <- move_peaks(mg_list_RT[["Ib_TRUE"]]
                                      , Sample = 'Fly_039'
                                      , peaks_list = movements_Ib_TRUE39 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ib_TRUE39 %>%
                                        pull(movement_dirs))
mg_list_area[["Ib_TRUE"]] <- move_peaks(mg_list_area[["Ib_TRUE"]]
                                        , Sample = 'Fly_039'
                                        , peaks_list = movements_Ib_TRUE39 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ib_TRUE39 %>%
                                          pull(movement_dirs))
# Sample 40
mg_list_RT[["Ib_TRUE"]] <- move_peaks(mg_list_RT[["Ib_TRUE"]]
                                      , Sample = 'Fly_040'
                                      , peaks_list = movements_Ib_TRUE40 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ib_TRUE40 %>%
                                        pull(movement_dirs))
mg_list_area[["Ib_TRUE"]] <- move_peaks(mg_list_area[["Ib_TRUE"]]
                                        , Sample = 'Fly_040'
                                        , peaks_list = movements_Ib_TRUE40 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ib_TRUE40 %>%
                                          pull(movement_dirs))
# Sample 33
mg_list_RT[["Ib_TRUE"]] <- move_peaks(mg_list_RT[["Ib_TRUE"]]
                                      , Sample = 'Fly_033'
                                      , peaks_list = movements_Ib_TRUE33 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ib_TRUE33 %>%
                                        pull(movement_dirs))
mg_list_area[["Ib_TRUE"]] <- move_peaks(mg_list_area[["Ib_TRUE"]]
                                        , Sample = 'Fly_033'
                                        , peaks_list = movements_Ib_TRUE33 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ib_TRUE33 %>%
                                          pull(movement_dirs))

# Sample 37
mg_list_RT[["Ib_TRUE"]] <- move_peaks(mg_list_RT[["Ib_TRUE"]]
                                      , Sample = 'Fly_037'
                                      , peaks_list = movements_Ib_TRUE37 %>%
                                        pull(peaks_list)
                                      , movement_dirs = movements_Ib_TRUE37 %>%
                                        pull(movement_dirs))
mg_list_area[["Ib_TRUE"]] <- move_peaks(mg_list_area[["Ib_TRUE"]]
                                        , Sample = 'Fly_037'
                                        , peaks_list = movements_Ib_TRUE37 %>%
                                          pull(peaks_list)
                                        , movement_dirs = movements_Ib_TRUE37 %>%
                                          pull(movement_dirs))

print("Alignment was corrected. Check heat map to verify that it is correct.")

# Check the corrected alignment ----
# This is done with the clustered heat map
## Data set with relative abundance of each peak (%)\
### Ca_FALSE
Ca_FALSE_per <- mg_list_area[['Ca_FALSE']] / rowSums(mg_list_area[['Ca_FALSE']]) * 100

### Ca_TRUE
Ca_TRUE_per <- mg_list_area[['Ca_TRUE']] / rowSums(mg_list_area[['Ca_TRUE']]) * 100

### Ib_FALSE
Ib_FALSE_per <- mg_list_area[['Ib_FALSE']] / rowSums(mg_list_area[['Ib_FALSE']]) * 100

### Ib_TRUE
Ib_TRUE_per <- mg_list_area[['Ib_TRUE']] / rowSums(mg_list_area[['Ib_TRUE']]) * 100

#### Set color palette
heatmap_colors <- viridis::turbo(200)

# set.seed(12345)
pdf(here("output", analysis, 'corrected-alignment_heatmap.pdf')
    , width = 20, height = 10)
{gplots::heatmap.2(as.matrix(Ca_FALSE_per %>% log1p()), 
                   #Rowv = F,
                   main = "Clustering of Carnica non-flying honeybee workers CHC",
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

{gplots::heatmap.2(as.matrix(Ca_TRUE_per %>% log1p()), 
                   main = "Clustering of Carnica flying honeybee workers CHC",
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

{gplots::heatmap.2(as.matrix(Ib_FALSE_per %>% log1p()), 
                   #Rowv = F,
                   main = "Clustering of Iberiensis non-flying honeybee workers CHC",
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

{gplots::heatmap.2(as.matrix(Ib_TRUE_per %>% log1p()), 
                   main = "Clustering of Iberiensis flying honeybee workers CHC",
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
### Ca_FALSE
mg_list_RT[["Ca_FALSE"]] <- mg_list_RT[["Ca_FALSE"]] %>% 
  select(-all_of(colnames(mg_list_RT[["Ca_FALSE"]][colSums(mg_list_RT[["Ca_FALSE"]]) == 0])))

mg_list_area[["Ca_FALSE"]] <- mg_list_area[["Ca_FALSE"]] %>% 
  select(-all_of(colnames(mg_list_area[["Ca_FALSE"]][colSums(mg_list_area[["Ca_FALSE"]]) == 0])))


mg_list_mean_RT[["Ca_FALSE"]] <- mg_list_mean_RT[["Ca_FALSE"]] %>% t %>% as.data.frame %>% 
  select(all_of(colnames(mg_list_RT[["Ca_FALSE"]]))) %>% t %>% as.data.frame

### Ca_TRUE
mg_list_RT[["Ca_TRUE"]] <- mg_list_RT[["Ca_TRUE"]] %>% 
  select(-all_of(colnames(mg_list_RT[["Ca_TRUE"]][colSums(mg_list_RT[["Ca_TRUE"]]) == 0])))

mg_list_area[["Ca_TRUE"]] <- mg_list_area[["Ca_TRUE"]] %>% 
  select(-all_of(colnames(mg_list_area[["Ca_TRUE"]][colSums(mg_list_area[["Ca_TRUE"]]) == 0])))


mg_list_mean_RT[["Ca_TRUE"]] <- mg_list_mean_RT[["Ca_TRUE"]] %>% t %>% as.data.frame %>% 
  select(all_of(colnames(mg_list_RT[["Ca_TRUE"]]))) %>% t %>% as.data.frame

### Ib_FALSE
mg_list_RT[["Ib_FALSE"]] <- mg_list_RT[["Ib_FALSE"]] %>% 
  select(-all_of(colnames(mg_list_RT[["Ib_FALSE"]][colSums(mg_list_RT[["Ib_FALSE"]]) == 0])))

mg_list_area[["Ib_FALSE"]] <- mg_list_area[["Ib_FALSE"]] %>% 
  select(-all_of(colnames(mg_list_area[["Ib_FALSE"]][colSums(mg_list_area[["Ib_FALSE"]]) == 0])))


mg_list_mean_RT[["Ib_FALSE"]] <- mg_list_mean_RT[["Ib_FALSE"]] %>% t %>% as.data.frame %>% 
  select(all_of(colnames(mg_list_RT[["Ib_FALSE"]]))) %>% t %>% as.data.frame

### Ib_TRUE
mg_list_RT[["Ib_TRUE"]] <- mg_list_RT[["Ib_TRUE"]] %>% 
  select(-all_of(colnames(mg_list_RT[["Ib_TRUE"]][colSums(mg_list_RT[["Ib_TRUE"]]) == 0])))

mg_list_area[["Ib_TRUE"]] <- mg_list_area[["Ib_TRUE"]] %>% 
  select(-all_of(colnames(mg_list_area[["Ib_TRUE"]][colSums(mg_list_area[["Ib_TRUE"]]) == 0])))


mg_list_mean_RT[["Ib_TRUE"]] <- mg_list_mean_RT[["Ib_TRUE"]] %>% t %>% as.data.frame %>% 
  select(all_of(colnames(mg_list_RT[["Ib_TRUE"]]))) %>% t %>% as.data.frame


### STD
STD_RT <- STD_RT %>% 
  select(-all_of(colnames(STD_RT[colSums(STD_RT) == 0])))

STD_mean_RT <- STD_mean_RT %>% t %>% as.data.frame %>% 
  select(all_of(colnames(STD_RT))) %>% t %>% as.data.frame

## Recalculate the mean RT
### Ca_FALSE
mg_list_RT[["Ca_FALSE"]][mg_list_RT[["Ca_FALSE"]] == 0] <- NA
mg_list_mean_RT[["Ca_FALSE"]] <- mg_list_mean_RT[["Ca_FALSE"]] %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = mg_list_RT[["Ca_FALSE"]] %>% colMeans(na.rm = T))

### Ca_TRUE
mg_list_RT[["Ca_TRUE"]][mg_list_RT[["Ca_TRUE"]] == 0] <- NA
mg_list_mean_RT[["Ca_TRUE"]] <- mg_list_mean_RT[["Ca_TRUE"]] %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = mg_list_RT[["Ca_TRUE"]] %>% colMeans(na.rm = T))

### Ib_FALSE
mg_list_RT[["Ib_FALSE"]][mg_list_RT[["Ib_FALSE"]] == 0] <- NA
mg_list_mean_RT[["Ib_FALSE"]] <- mg_list_mean_RT[["Ib_FALSE"]] %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = mg_list_RT[["Ib_FALSE"]] %>% colMeans(na.rm = T))

### Ib_TRUE
mg_list_RT[["Ib_TRUE"]][mg_list_RT[["Ib_TRUE"]] == 0] <- NA
mg_list_mean_RT[["Ib_TRUE"]] <- mg_list_mean_RT[["Ib_TRUE"]] %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = mg_list_RT[["Ib_TRUE"]] %>% colMeans(na.rm = T))

### STD
STD_RT[STD_RT == 0] <- NA
STD_mean_RT <- STD_mean_RT %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = STD_RT %>% colMeans(na.rm = T))


## re-shape data frames before exporting them
### Ca_FALSE
mg_list_RT[['Ca_FALSE']] <- cbind.data.frame(Peak = row.names(mg_list_mean_RT[['Ca_FALSE']])
                               , mean_RT = mg_list_mean_RT[['Ca_FALSE']]$new_mean_RT
                               # Place Ca_FALSE as columns
                               , mg_list_RT[['Ca_FALSE']] %>%
                                 t %>% 
                                 as.data.frame %>% 
                                 # Order them in descending order
                                 # , regarding their total abundance
                                 select(all_of(mg_list_area[['Ca_FALSE']] %>% 
                                                 rowSums() %>% 
                                                 sort(decreasing = T) %>%
                                                 names()))) %>%
  as_tibble()
mg_list_RT[['Ca_FALSE']]

mg_list_area[['Ca_FALSE']] <- cbind.data.frame(Peak = row.names(mg_list_mean_RT[['Ca_FALSE']])
                                 , mean_RT = mg_list_mean_RT[['Ca_FALSE']]$new_mean_RT
                                 # Place Ca_FALSE as columns
                                 , mg_list_area[['Ca_FALSE']] %>%
                                   t %>% 
                                   as.data.frame %>% 
                                   # Order them in descending order
                                   # , regarding their total abundance
                                   select(all_of(mg_list_area[['Ca_FALSE']] %>% 
                                                   rowSums() %>% 
                                                   sort(decreasing = T) %>%
                                                   names()))) %>%
  as_tibble()
mg_list_area[['Ca_FALSE']]

### Ca_TRUE
mg_list_RT[['Ca_TRUE']] <- cbind.data.frame(Peak = row.names(mg_list_mean_RT[['Ca_TRUE']])
                                             , mean_RT = mg_list_mean_RT[['Ca_TRUE']]$new_mean_RT
                                             # Place Ca_TRUE as columns
                                             , mg_list_RT[['Ca_TRUE']] %>%
                                               t %>% 
                                               as.data.frame %>% 
                                               # Order them in descending order
                                               # , regarding their total abundance
                                               select(all_of(mg_list_area[['Ca_TRUE']] %>% 
                                                               rowSums() %>% 
                                                               sort(decreasing = T) %>%
                                                               names()))) %>%
  as_tibble()
mg_list_RT[['Ca_TRUE']]

mg_list_area[['Ca_TRUE']] <- cbind.data.frame(Peak = row.names(mg_list_mean_RT[['Ca_TRUE']])
                                               , mean_RT = mg_list_mean_RT[['Ca_TRUE']]$new_mean_RT
                                               # Place Ca_TRUE as columns
                                               , mg_list_area[['Ca_TRUE']] %>%
                                                 t %>% 
                                                 as.data.frame %>% 
                                                 # Order them in descending order
                                                 # , regarding their total abundance
                                                 select(all_of(mg_list_area[['Ca_TRUE']] %>% 
                                                                 rowSums() %>% 
                                                                 sort(decreasing = T) %>%
                                                                 names()))) %>%
  as_tibble()
mg_list_area[['Ca_TRUE']]

### Ib_FALSE
mg_list_RT[['Ib_FALSE']] <- cbind.data.frame(Peak = row.names(mg_list_mean_RT[['Ib_FALSE']])
                                             , mean_RT = mg_list_mean_RT[['Ib_FALSE']]$new_mean_RT
                                             # Place Ib_FALSE as columns
                                             , mg_list_RT[['Ib_FALSE']] %>%
                                               t %>% 
                                               as.data.frame %>% 
                                               # Order them in descending order
                                               # , regarding their total abundance
                                               select(all_of(mg_list_area[['Ib_FALSE']] %>% 
                                                               rowSums() %>% 
                                                               sort(decreasing = T) %>%
                                                               names()))) %>%
  as_tibble()
mg_list_RT[['Ib_FALSE']]

mg_list_area[['Ib_FALSE']] <- cbind.data.frame(Peak = row.names(mg_list_mean_RT[['Ib_FALSE']])
                                               , mean_RT = mg_list_mean_RT[['Ib_FALSE']]$new_mean_RT
                                               # Place Ib_FALSE as columns
                                               , mg_list_area[['Ib_FALSE']] %>%
                                                 t %>% 
                                                 as.data.frame %>% 
                                                 # Order them in descending order
                                                 # , regarding their total abundance
                                                 select(all_of(mg_list_area[['Ib_FALSE']] %>% 
                                                                 rowSums() %>% 
                                                                 sort(decreasing = T) %>%
                                                                 names()))) %>%
  as_tibble()
mg_list_area[['Ib_FALSE']]

### Ib_TRUE
mg_list_RT[['Ib_TRUE']] <- cbind.data.frame(Peak = row.names(mg_list_mean_RT[['Ib_TRUE']])
                                            , mean_RT = mg_list_mean_RT[['Ib_TRUE']]$new_mean_RT
                                            # Place Ib_TRUE as columns
                                            , mg_list_RT[['Ib_TRUE']] %>%
                                              t %>% 
                                              as.data.frame %>% 
                                              # Order them in descending order
                                              # , regarding their total abundance
                                              select(all_of(mg_list_area[['Ib_TRUE']] %>% 
                                                              rowSums() %>% 
                                                              sort(decreasing = T) %>%
                                                              names()))) %>%
  as_tibble()
mg_list_RT[['Ib_TRUE']]

mg_list_area[['Ib_TRUE']] <- cbind.data.frame(Peak = row.names(mg_list_mean_RT[['Ib_TRUE']])
                                              , mean_RT = mg_list_mean_RT[['Ib_TRUE']]$new_mean_RT
                                              # Place Ib_TRUE as columns
                                              , mg_list_area[['Ib_TRUE']] %>%
                                                t %>% 
                                                as.data.frame %>% 
                                                # Order them in descending order
                                                # , regarding their total abundance
                                                select(all_of(mg_list_area[['Ib_TRUE']] %>% 
                                                                rowSums() %>% 
                                                                sort(decreasing = T) %>%
                                                                names()))) %>%
  as_tibble()
mg_list_area[['Ib_TRUE']]

### STD
STD_RT <- cbind.data.frame(Peak = row.names(STD_mean_RT)
                             , mean_RT = STD_mean_RT$new_mean_RT
                             # Place STD as columns
                             , STD_RT %>%
                               t %>% 
                               as.data.frame) %>% 
  as_tibble()
STD_RT



## Save the data frames 
save(list = c("mg_list_area", "mg_list_RT"
            , "STD_RT")
     , file = here("data"
                   , "processed", analysis, "aligned_gcms-data.Rdata"))
print("The aligned data frames were exported")

DF_names <- names(mg_list_RT) 
for(df in DF_names){
  name_export <- df
  openxlsx::write.xlsx(mg_list_RT[[df]]
                       , here("data"
                              , "raw"
                              , analysis
                              , "tmp"
                              , paste0(name_export
                                              ,"_RT_"
                                              ,"table.xlsx")))
}

openxlsx::write.xlsx(STD_RT
                     , here("data", "raw", analysis
                            , "tmp", "STD_RT_table.xlsx"))

print("The RT tables for CHC identification was exported")

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", analysis, "SInf_Script03.txt"))
print("The sessionInfo report was exported. The script 03 finished running")

# ## Detach/unload packages
# lapply(NPacks, unloadNamespace)
# 
# ## Clear environment
# rm(list=ls())


