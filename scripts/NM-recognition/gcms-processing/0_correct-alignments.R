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
STD_0322_RT <- STD_0322_RT %>% 
  select(-mean_RT) %>% 
  t %>% 
  as.data.frame()
str(STD_0322_RT)

### STD_0422
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
                                                     , 187, 189)))
                               , movement_dirs = c('up', 'up', 'up', 'up'
                                                 , 'up', 'up', 'up', 'up'
                                                 , 'up'))
str(movements_STD_0422_H2204_2)







## Move peaks within samples to correct alignments ####
### The function move_peaks() displaces the specified peaks of a sample 
### by one position up (previous table peak)/down (later table peak)
### It has to be applied to both RT and area data frames
### Br_Ca ####
# Sample 194
Br_Ca_RT <- move_peaks(Br_Ca_RT
                         , Sample = 'NM_194'
                         , peaks_list = movements_Br_Ca194 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Br_Ca194 %>%
                           pull(movement_dirs))
Br_Ca_area <- move_peaks(Br_Ca_area
                           , Sample = 'NM_194'
                           , peaks_list = movements_Br_Ca194 %>%
                             pull(peaks_list)
                           , movement_dirs = movements_Br_Ca194 %>%
                             pull(movement_dirs))

# Sample 197
Br_Ca_RT <- move_peaks(Br_Ca_RT
                       , Sample = 'NM_197'
                       , peaks_list = movements_Br_Ca197 %>%
                         pull(peaks_list)
                       , movement_dirs = movements_Br_Ca197 %>%
                         pull(movement_dirs))
Br_Ca_area <- move_peaks(Br_Ca_area
                         , Sample = 'NM_197'
                         , peaks_list = movements_Br_Ca197 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Br_Ca197 %>%
                           pull(movement_dirs))

### Wu_Ca ####
# Sample 63
Wu_Ca_RT <- move_peaks(Wu_Ca_RT
                       , Sample = 'NM_063_1'
                       , peaks_list = movements_Wu_Ca63 %>%
                         pull(peaks_list)
                       , movement_dirs = movements_Wu_Ca63 %>%
                         pull(movement_dirs))
Wu_Ca_area <- move_peaks(Wu_Ca_area
                         , Sample = 'NM_063_1'
                         , peaks_list = movements_Wu_Ca63 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Wu_Ca63 %>%
                           pull(movement_dirs))

# Sample 64
Wu_Ca_RT <- move_peaks(Wu_Ca_RT
                       , Sample = 'NM_064_1'
                       , peaks_list = movements_Wu_Ca64 %>%
                         pull(peaks_list)
                       , movement_dirs = movements_Wu_Ca64 %>%
                         pull(movement_dirs))
Wu_Ca_area <- move_peaks(Wu_Ca_area
                         , Sample = 'NM_064_1'
                         , peaks_list = movements_Wu_Ca64 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Wu_Ca64 %>%
                           pull(movement_dirs))

# Sample 87
Wu_Ca_RT <- move_peaks(Wu_Ca_RT
                       , Sample = 'NM_087'
                       , peaks_list = movements_Wu_Ca87 %>%
                         pull(peaks_list)
                       , movement_dirs = movements_Wu_Ca87 %>%
                         pull(movement_dirs))
Wu_Ca_area <- move_peaks(Wu_Ca_area
                         , Sample = 'NM_087'
                         , peaks_list = movements_Wu_Ca87 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Wu_Ca87 %>%
                           pull(movement_dirs))

# Sample 94
Wu_Ca_RT <- move_peaks(Wu_Ca_RT
                       , Sample = 'NM_094'
                       , peaks_list = movements_Wu_Ca94 %>%
                         pull(peaks_list)
                       , movement_dirs = movements_Wu_Ca94 %>%
                         pull(movement_dirs))
Wu_Ca_area <- move_peaks(Wu_Ca_area
                         , Sample = 'NM_094'
                         , peaks_list = movements_Wu_Ca94 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Wu_Ca94 %>%
                           pull(movement_dirs))

# Sample 166
Wu_Ca_RT <- move_peaks(Wu_Ca_RT
                       , Sample = 'NM_166'
                       , peaks_list = movements_Wu_Ca166 %>%
                         pull(peaks_list)
                       , movement_dirs = movements_Wu_Ca166 %>%
                         pull(movement_dirs))
Wu_Ca_area <- move_peaks(Wu_Ca_area
                         , Sample = 'NM_166'
                         , peaks_list = movements_Wu_Ca166 %>%
                           pull(peaks_list)
                         , movement_dirs = movements_Wu_Ca166 %>%
                           pull(movement_dirs))
### STD_0422 ####
# Sample L_SSL_2204_2
STD_0422_RT <- move_peaks(STD_0422_RT
                       , Sample = 'L_SSL_2204_2'
                       , peaks_list = movements_STD_0422_L2204_2 %>%
                         pull(peaks_list)
                       , movement_dirs = movements_STD_0422_L2204_2 %>%
                         pull(movement_dirs))

# Sample H_SSL_2204_2
STD_0422_RT <- move_peaks(STD_0422_RT
                          , Sample = 'H_SSL_2204_2'
                          , peaks_list = movements_STD_0422_H2204_2 %>%
                            pull(peaks_list)
                          , movement_dirs = movements_STD_0422_H2204_2 %>%
                            pull(movement_dirs))

print("Alignment was corrected. Check heat map to verify that it is correct.")

# Check the corrected alignment ----
# This is done with the clustered heat map
## Data set with relative abundance of each peak (%)\
### Br_Ca
Br_Ca_per <- Br_Ca_area / rowSums(Br_Ca_area) * 100

### Br_Ib
Br_Ib_per <- Br_Ib_area / rowSums(Br_Ib_area) * 100

### Wu_Ca
Wu_Ca_per <- Wu_Ca_area / rowSums(Wu_Ca_area) * 100

### Wu_Ib_0322
Wu_Ib_0322_per <- Wu_Ib_0322_area / rowSums(Wu_Ib_0322_area) * 100

### Wu_Ib_0422
Wu_Ib_0422_per <- Wu_Ib_0422_area / rowSums(Wu_Ib_0422_area) * 100

#### Set color palette
heatmap_colors <- viridis::turbo(200)

# set.seed(12345)
pdf(here("output", "NM-Recognition", 'corrected-alginment_heatmap.pdf')
    , width = 20, height = 10)
{gplots::heatmap.2(as.matrix(Br_Ca_per %>% log1p()), 
                   #Rowv = F,
                   main = "Clustering of Braganza Carnica honeybee workers CHC",
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

{gplots::heatmap.2(as.matrix(Br_Ib_per %>% log1p()), 
                   main = "Clustering of Braganza Iberiensis honeybee workers 
                   CHC",
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

{gplots::heatmap.2(as.matrix(Wu_Ca_per %>% log1p()), 
                   main = "Clustering of Wuerzburg Carnica honeybee workers 
                   CHC",
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

{gplots::heatmap.2(as.matrix(Wu_Ib_0322_per %>% log1p()), 
                   main = "Clustering of Wuerzburg Iberiensis (03/22) honeybee 
                   workers CHC",
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

{gplots::heatmap.2(as.matrix(Wu_Ib_0422_per %>% log1p()), 
                   main = "Clustering of Wuerzburg Iberiensis (04/22) honeybee 
                   workers CHC",
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
# 
# 
# {gplots::heatmap.2(as.matrix(STD_0322_per %>% log1p()), 
#                    main = "Clustering of standard alkanes from 03/22", 
#                    srtCol = 90,
#                    dendrogram = "row",
#                    # Rowv = dend,
#                    Colv = "NA", # this to make sure the columns are not ordered
#                    trace = "none",          
#                    # margins =c(5,0.1),      
#                    key.xlab = "Relative abundance (%) on Cuticle",
#                    #denscol = "grey",
#                    #density.info = "density",
#                    # RowSideColors = microlab, # to add nice colored strips        
#                    col = heatmap_colors)}
# 
# {gplots::heatmap.2(as.matrix(STD_0422_per %>% log1p()), 
#                    main = "Clustering of standard alkanes from 04/22", 
#                    srtCol = 90,
#                    dendrogram = "row",
#                    # Rowv = dend,
#                    Colv = "NA", # this to make sure the columns are not ordered
#                    trace = "none",          
#                    # margins =c(5,0.1),      
#                    key.xlab = "Relative abundance (%) on Cuticle",
#                    #denscol = "grey",
#                    #density.info = "density",
#                    # RowSideColors = microlab, # to add nice colored strips        
#                    col = heatmap_colors)}
dev.off()
print("Heat map to verify corrected alignment was exported")

# Export the aligned data for CHC identification ----
## Drop empty columns
### Br_Ca
Br_Ca_RT <- Br_Ca_RT %>% 
  select(-all_of(colnames(Br_Ca_RT[colSums(Br_Ca_RT) == 0])))

Br_Ca_area <- Br_Ca_area %>% 
  select(-all_of(colnames(Br_Ca_area[colSums(Br_Ca_area) == 0])))

Br_Ca_mean_RT <- Br_Ca_mean_RT %>% t %>% as.data.frame %>% 
  select(all_of(colnames(Br_Ca_RT))) %>% t %>% as.data.frame

### Br_Ib
Br_Ib_RT <- Br_Ib_RT %>% 
  select(-all_of(colnames(Br_Ib_RT[colSums(Br_Ib_RT) == 0])))

Br_Ib_area <- Br_Ib_area %>% 
  select(-all_of(colnames(Br_Ib_area[colSums(Br_Ib_area) == 0])))

Br_Ib_mean_RT <- Br_Ib_mean_RT %>% t %>% as.data.frame %>% 
  select(all_of(colnames(Br_Ib_RT))) %>% t %>% as.data.frame

### Wu_Ca
Wu_Ca_RT <- Wu_Ca_RT %>% 
  select(-all_of(colnames(Wu_Ca_RT[colSums(Wu_Ca_RT) == 0])))

Wu_Ca_area <- Wu_Ca_area %>% 
  select(-all_of(colnames(Wu_Ca_area[colSums(Wu_Ca_area) == 0])))

Wu_Ca_mean_RT <- Wu_Ca_mean_RT %>% t %>% as.data.frame %>% 
  select(all_of(colnames(Wu_Ca_RT))) %>% t %>% as.data.frame


### Wu_Ib_0322
Wu_Ib_0322_RT <- Wu_Ib_0322_RT %>% 
  select(-all_of(colnames(Wu_Ib_0322_RT[colSums(Wu_Ib_0322_RT) == 0])))

Wu_Ib_0322_area <- Wu_Ib_0322_area %>% 
  select(-all_of(colnames(Wu_Ib_0322_area[colSums(Wu_Ib_0322_area) == 0])))

Wu_Ib_0322_mean_RT <- Wu_Ib_0322_mean_RT %>% t %>% as.data.frame %>% 
  select(all_of(colnames(Wu_Ib_0322_RT))) %>% t %>% as.data.frame

### Wu_Ib_0422
Wu_Ib_0422_RT <- Wu_Ib_0422_RT %>% 
  select(-all_of(colnames(Wu_Ib_0422_RT[colSums(Wu_Ib_0422_RT) == 0])))

Wu_Ib_0422_area <- Wu_Ib_0422_area %>% 
  select(-all_of(colnames(Wu_Ib_0422_area[colSums(Wu_Ib_0422_area) == 0])))

Wu_Ib_0422_mean_RT <- Wu_Ib_0422_mean_RT %>% t %>% as.data.frame %>% 
  select(all_of(colnames(Wu_Ib_0422_RT))) %>% t %>% as.data.frame


### STD_0322
STD_0322_RT <- STD_0322_RT %>% 
  select(-all_of(colnames(STD_0322_RT[colSums(STD_0322_RT) == 0])))

STD_0322_mean_RT <- STD_0322_mean_RT %>% t %>% as.data.frame %>% 
  select(all_of(colnames(STD_0322_RT))) %>% t %>% as.data.frame

### STD_0422
STD_0422_RT <- STD_0422_RT %>% 
  select(-all_of(colnames(STD_0422_RT[colSums(STD_0422_RT) == 0])))

STD_0422_mean_RT <- STD_0422_mean_RT %>% t %>% as.data.frame %>% 
  select(all_of(colnames(STD_0422_RT))) %>% t %>% as.data.frame


## Recalculate the mean RT
### Br_Ca
Br_Ca_RT[Br_Ca_RT == 0] <- NA
Br_Ca_mean_RT <- Br_Ca_mean_RT %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = Br_Ca_RT %>% colMeans(na.rm = T))

### Br_Ib
Br_Ib_RT[Br_Ib_RT == 0] <- NA
Br_Ib_mean_RT <- Br_Ib_mean_RT %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = Br_Ib_RT %>% colMeans(na.rm = T))

### Wu_Ca
Wu_Ca_RT[Wu_Ca_RT == 0] <- NA
Wu_Ca_mean_RT <- Wu_Ca_mean_RT %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = Wu_Ca_RT %>% colMeans(na.rm = T))

### Wu_Ib_0322
Wu_Ib_0322_RT[Wu_Ib_0322_RT == 0] <- NA
Wu_Ib_0322_mean_RT <- Wu_Ib_0322_mean_RT %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = Wu_Ib_0322_RT %>% colMeans(na.rm = T))

### Wu_Ib_0422
Wu_Ib_0422_RT[Wu_Ib_0422_RT == 0] <- NA
Wu_Ib_0422_mean_RT <- Wu_Ib_0422_mean_RT %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = Wu_Ib_0422_RT %>% colMeans(na.rm = T))

### STD_0422
STD_0422_RT[STD_0422_RT == 0] <- NA
STD_0422_mean_RT <- STD_0422_mean_RT %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = STD_0422_RT %>% colMeans(na.rm = T))


### STD_0322
STD_0322_RT[STD_0322_RT == 0] <- NA
STD_0322_mean_RT <- STD_0322_mean_RT %>% 
  transmute(old_mean_RT = mean_RT) %>% 
  mutate(new_mean_RT = STD_0322_RT %>% colMeans(na.rm = T))

## re-shape data frames before exporting them
### Br_Ca
Br_Ca_RT <- cbind.data.frame(Peak = row.names(Br_Ca_mean_RT)
                               , mean_RT = Br_Ca_mean_RT$new_mean_RT
                               # Place Br_Ca as columns
                               , Br_Ca_RT %>%
                                 t %>% 
                                 as.data.frame %>% 
                                 # Order them in descending order
                                 # , regarding their total abundance
                                 select(all_of(Br_Ca_area %>% 
                                                 rowSums() %>% 
                                                 sort(decreasing = T) %>%
                                                 names()))) %>%
  as_tibble()
Br_Ca_RT

Br_Ca_area <- cbind.data.frame(Peak = row.names(Br_Ca_mean_RT)
                                 , mean_RT = Br_Ca_mean_RT$new_mean_RT
                                 # Place Br_Ca as columns
                                 , Br_Ca_area %>%
                                   t %>% 
                                   as.data.frame %>% 
                                   # Order them in descending order
                                   # , regarding their total abundance
                                   select(all_of(Br_Ca_area %>% 
                                                   rowSums() %>% 
                                                   sort(decreasing = T) %>%
                                                   names()))) %>%
  as_tibble()
Br_Ca_area

### Br_Ib
Br_Ib_RT <- cbind.data.frame(Peak = row.names(Br_Ib_mean_RT)
                             , mean_RT = Br_Ib_mean_RT$new_mean_RT
                             # Place Br_Ib as columns
                             , Br_Ib_RT %>%
                               t %>% 
                               as.data.frame %>% 
                               # Order them in descending order
                               # , regarding their total abundance
                               select(all_of(Br_Ib_area %>% 
                                               rowSums() %>% 
                                               sort(decreasing = T) %>%
                                               names()))) %>%
  as_tibble()
Br_Ib_RT

Br_Ib_area <- cbind.data.frame(Peak = row.names(Br_Ib_mean_RT)
                               , mean_RT = Br_Ib_mean_RT$new_mean_RT
                               # Place Br_Ib as columns
                               , Br_Ib_area %>%
                                 t %>% 
                                 as.data.frame %>% 
                                 # Order them in descending order
                                 # , regarding their total abundance
                                 select(all_of(Br_Ib_area %>% 
                                                 rowSums() %>% 
                                                 sort(decreasing = T) %>%
                                                 names()))) %>%
  as_tibble()
Br_Ib_area

### Wu_Ca
Wu_Ca_RT <- cbind.data.frame(Peak = row.names(Wu_Ca_mean_RT)
                             , mean_RT = Wu_Ca_mean_RT$new_mean_RT
                             # Place Wu_Ca as columns
                             , Wu_Ca_RT %>%
                               t %>% 
                               as.data.frame %>% 
                               # Order them in descending order
                               # , regarding their total abundance
                               select(all_of(Wu_Ca_area %>% 
                                               rowSums() %>% 
                                               sort(decreasing = T) %>%
                                               names()))) %>%
  as_tibble()
Wu_Ca_RT

Wu_Ca_area <- cbind.data.frame(Peak = row.names(Wu_Ca_mean_RT)
                               , mean_RT = Wu_Ca_mean_RT$new_mean_RT
                               # Place Wu_Ca as columns
                               , Wu_Ca_area %>%
                                 t %>% 
                                 as.data.frame %>% 
                                 # Order them in descending order
                                 # , regarding their total abundance
                                 select(all_of(Wu_Ca_area %>% 
                                                 rowSums() %>% 
                                                 sort(decreasing = T) %>%
                                                 names()))) %>%
  as_tibble()
Wu_Ca_area

### Wu_Ib_0322
Wu_Ib_0322_RT <- cbind.data.frame(Peak = row.names(Wu_Ib_0322_mean_RT)
                             , mean_RT = Wu_Ib_0322_mean_RT$new_mean_RT
                             # Place Wu_Ib_0322 as columns
                             , Wu_Ib_0322_RT %>%
                               t %>% 
                               as.data.frame %>% 
                               # Order them in descending order
                               # , regarding their total abundance
                               select(all_of(Wu_Ib_0322_area %>% 
                                               rowSums() %>% 
                                               sort(decreasing = T) %>%
                                               names()))) %>%
  as_tibble()
Wu_Ib_0322_RT

Wu_Ib_0322_area <- cbind.data.frame(Peak = row.names(Wu_Ib_0322_mean_RT)
                               , mean_RT = Wu_Ib_0322_mean_RT$new_mean_RT
                               # Place Wu_Ib_0322 as columns
                               , Wu_Ib_0322_area %>%
                                 t %>% 
                                 as.data.frame %>% 
                                 # Order them in descending order
                                 # , regarding their total abundance
                                 select(all_of(Wu_Ib_0322_area %>% 
                                                 rowSums() %>% 
                                                 sort(decreasing = T) %>%
                                                 names()))) %>%
  as_tibble()
Wu_Ib_0322_area

### Wu_Ib_0422
Wu_Ib_0422_RT <- cbind.data.frame(Peak = row.names(Wu_Ib_0422_mean_RT)
                             , mean_RT = Wu_Ib_0422_mean_RT$new_mean_RT
                             # Place Wu_Ib_0422 as columns
                             , Wu_Ib_0422_RT %>%
                               t %>% 
                               as.data.frame %>% 
                               # Order them in descending order
                               # , regarding their total abundance
                               select(all_of(Wu_Ib_0422_area %>% 
                                               rowSums() %>% 
                                               sort(decreasing = T) %>%
                                               names()))) %>%
  as_tibble()
Wu_Ib_0422_RT

Wu_Ib_0422_area <- cbind.data.frame(Peak = row.names(Wu_Ib_0422_mean_RT)
                               , mean_RT = Wu_Ib_0422_mean_RT$new_mean_RT
                               # Place Wu_Ib_0422 as columns
                               , Wu_Ib_0422_area %>%
                                 t %>% 
                                 as.data.frame %>% 
                                 # Order them in descending order
                                 # , regarding their total abundance
                                 select(all_of(Wu_Ib_0422_area %>% 
                                                 rowSums() %>% 
                                                 sort(decreasing = T) %>%
                                                 names()))) %>%
  as_tibble()
Wu_Ib_0422_area

### STD_0422
STD_0422_RT <- cbind.data.frame(Peak = row.names(STD_0422_mean_RT)
                             , mean_RT = STD_0422_mean_RT$new_mean_RT
                             # Place STD_0422 as columns
                             , STD_0422_RT %>%
                               t %>% 
                               as.data.frame) %>% 
  as_tibble()
STD_0422_RT


### STD_0322
STD_0322_RT <- cbind.data.frame(Peak = row.names(STD_0322_mean_RT)
                             , mean_RT = STD_0322_mean_RT$new_mean_RT
                             # Place STD_0322 as columns
                             , STD_0322_RT %>%
                               t %>% 
                               as.data.frame) %>%
  as_tibble()
STD_0322_RT



## Save the data frames 
save(list = c("Br_Ca_area", "Br_Ca_RT", "Br_Ib_area", "Br_Ib_RT"
            , "Wu_Ca_area", "Wu_Ca_RT", "Wu_Ib_0322_area", "Wu_Ib_0322_RT"
            , "Wu_Ib_0422_area", "Wu_Ib_0422_RT"
            , "STD_0422_RT", "STD_0322_RT")
     , file = here("data"
                   , "processed", "NM-recognition", "aligned_gcms-data.Rdata"))
print("The aligned data frames were exported")

openxlsx::write.xlsx(Br_Ca_RT
                     , here("data", "raw", "NM-recognition"
                            , "tmp", "Br_Ca_RT_table.xlsx"))

openxlsx::write.xlsx(Br_Ib_RT
                     , here("data", "raw", "NM-recognition"
                            , "tmp", "Br_Ib_RT_table.xlsx"))

openxlsx::write.xlsx(Wu_Ca_RT
                     , here("data", "raw", "NM-recognition"
                            , "tmp", "Wu_Ca_RT_table.xlsx"))

openxlsx::write.xlsx(Wu_Ib_0322_RT
                     , here("data", "raw", "NM-recognition"
                            , "tmp", "Wu_Ib_0322_RT_table.xlsx"))

openxlsx::write.xlsx(Wu_Ib_0422_RT
                     , here("data", "raw", "NM-recognition"
                            , "tmp", "Wu_Ib_0422_RT_table.xlsx"))

openxlsx::write.xlsx(STD_0422_RT
                     , here("data", "raw", "NM-recognition"
                            , "tmp", "STD_0422_RT_table.xlsx"))

openxlsx::write.xlsx(STD_0322_RT
                     , here("data", "raw", "NM-recognition"
                            , "tmp", "STD_0322_RT_table.xlsx"))


print("The RT tables for CHC identification was exported")

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "NM-recognition", "SInf_Script03.txt"))
print("The sessionInfo report was exported. The script 03 finished running")

# ## Detach/unload packages
# lapply(NPacks, unloadNamespace)
# 
# ## Clear environment
# rm(list=ls())


