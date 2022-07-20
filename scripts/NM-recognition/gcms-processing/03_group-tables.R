# Install pacman package if necessary ----
#install.packages("pacman")

# Load required packages ----
NPacks <- c("readr", "stringr","tidyr", "dplyr", "here")
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

## Function to remove trace compounds from each sample
trace_comps <- function(group_daten, threshold) {
  cat('\n')
  ### Verify the total abundance per sample before deleting trace compounds
  print("Total abundance per sample before deleting trace compounds")
  print(group_daten %>% colSums())
  
  ### Calculate the relative abundance (%) of each compound per sample 
  group_daten_percent <- group_daten %>% t / rowSums(group_daten %>% t) * 100
  group_daten_percent <- group_daten_percent %>% 
    t %>% 
    as.data.frame()
  
  cat('\n')
  ### Verify that the sum of all relative abundances per sample is exactly 100
  print("Total relative abundance per sample before deleting trace compounds")
  print(group_daten_percent %>% colSums())
  
  ### Delete every peak of a sample that represent less
  ### than the 0.01% of the sample
  group_daten[group_daten_percent < threshold] <- NA
  
  cat('\n')
  ### Verify the total abundance per sample after deleting trace compounds
  print("Total abundance per sample after deleting trace compounds")
  print(group_daten %>% colSums(na.rm = T))
  cat('\n')
  
  group_daten
}

## Function to remove, from each group data, the compounds that are rare to
## that group
below_50 <- function(df_area, grouping.info){
  require(dplyr)
  
  grouping.info <- grouping.info %>%
    unite('group', !Individual, remove = F)
  
  print(paste("Total number of groups"
              , length(unique(grouping.info$group))
              , sep = ": "))
  cat('\n')
  
  groups_list <- list()
  lap <- 1
  for (g in unique(grouping.info$group)) {
    print(paste(paste("Group No."
                      , lap
                      , sep = " ")
                , g
                , sep = ": "))
    print(paste("Samples of the group:"
                , grouping.info %>% 
                  filter(group == g) %>% 
                  pull(Individual) %>% 
                  as.character() %>% 
                  paste0(collapse = " ")
                , sep = " "))
    
    df_g <- df_area %>% 
      select(all_of(grouping.info %>% 
                      filter(group == g) %>% 
                      pull(Individual)))
    df_g <- cbind(df_area %>% 
                    select(Peak:Compound)
                  , df_g) %>% 
      as_tibble()
    
    # Loop to get the list of peaks containing compounds that are rare to the group
    rare_comps_g <- c()
    na_proportions <- c()
    for (p in df_g$Peak) {
      # Calculate the proportion of samples within the group that do not have the compound of
      # the corresponding peak, by calculating the number of NAs in the row of the peak
      # and dividing it by the number of samples in the group
      na_prop <- df_g %>% 
        filter(Peak == p) %>% 
        select(!Peak:Compound) %>% 
        is.na() %>% 
        sum() / df_g %>% 
        select(!Peak:Compound) %>% ncol()
      
      # Is the proportion of NAs higher than 50%? 
      # If yes, the peak contains a rare compound
      if(0.5 < na_prop){
        rare_comps_g <- c(rare_comps_g, p)
        na_proportions <- c(na_proportions, na_prop)
      }
    }
    
    cat('\n')
    print(paste("The following"
                , length(rare_comps_g)
                , "compounds are rare (frequency < 0.5) within"
                , g
                , "group"))
    print("Therefore, they will be deleted from the data set of that group")
    data.frame(compounds = rare_comps_g, frequency = 1 - na_proportions) %>% 
      print()
    cat('\n')
    
    # Erase, from the group data, the peaks containing compounds that are rare to the group  
    df_g <- df_g %>% filter(!Peak %in% rare_comps_g)
    
    # Store the modified group data frame into groups_list
    groups_list[[g]] <- df_g
    lap <- lap + 1
  }
  
  cat('\n')
  print(" All rare compounds were removed from the corresponding group data")
  cat('\n')
  
  # Loop to automatically merge all the data frames in groups_list into a single one
  new_df_area <- groups_list[[1]] %>% select(Peak:Compound)
  for (i in groups_list) {
    new_df_area <- my_merge(new_df_area, i, merge_by = c("Peak", "Compound"))
  }
  # Return the modified df_area
  new_df_area %>% as_tibble()
}

## Function to calculate the Kováts retention index of the compounds in
## the samples
k.ri <- function(Comps) {
  
  cat('\n')
  require(dplyr)
  
  Comps2 <- Comps # New data frame to avoid modifications in the original
  
  # list of peaks containing target compounds (identified CHC)
  comps.list <- Comps2 %>% 
    filter(Class != "STD") %>% 
    pull(Peak)
  
  print(paste("The RI of"
              , length(comps.list)
              , "compounds will be calculated"
              , sep = " "))
  
  # Standards are considered alkanes
  Comps2["Class"][Comps2["Class"] == "STD"] <- "Alkane"
  
  # A data frame listign the alkanes
  alka.list <- Comps2 %>%
    filter(Class == "Alkane") %>%
    pull(Compound)
  
  # Empty vector to store the retention indexes
  ri.list <- c()
  Lap = 1 # set the lap count on 1
  
  # Loop making iterations for each compound
  for (n in comps.list) {
    # Data frame with the info of the compound
    # that is defined by current iteration
    n.comp <- Comps2 %>% filter(Comps2$Peak == n)
    
    # Report lap number and compound
    cat('\n')
    print(paste("Compound", Lap, n.comp[, "Compound"], sep = " "))
    
    # Calculate the RI depending on whether the compound is or not an alkane
    if (!n.comp$Compound %in% alka.list) {
      
      # Report current compound name
      # and that it is not considered as an alkane
      print(paste(n.comp[, "Compound"], "is not an Alkane", sep = " "))
      
      # get rt of current compound
      rt <- n.comp %>% pull(mean_RT)
      
      # Extract previous alkane data
      prev.alka <- Comps2 %>%
        filter(mean_RT < rt) %>%
        filter(Class == "Alkane")
      prev.alka <- prev.alka %>%
        filter(Chain.length == max(prev.alka$Chain.length))
      
      # Extract next alkane data
      next.alka <- Comps2 %>%
        filter(mean_RT > rt) %>%
        filter(Class == "Alkane")
      next.alka <- next.alka %>%
        filter(Chain.length == min(next.alka$Chain.length))
      
      # Calculate the retention index for the target compound
      rit <- 100 *
        (((n.comp$mean_RT - prev.alka$mean_RT) /
            (next.alka$mean_RT - prev.alka$mean_RT)) +
           prev.alka$Chain.length)
      
      # Round up the RI to make it an integer
      rit <- rit %>% round(digits = 0)
      ri.list <- c(ri.list, rit) # Store the RI in ri.list
      print(paste("RI:", rit, sep = " ")) # Report RI
    } else {
      # Report current compound/peak name
      # and that it is considered as an alkane
      print(paste(n.comp[, "Compound"], "is an Alkane", sep = " "))
      
      # Calculate the retention index for the target compound
      rit <- n.comp$Chain.length * 100
      ri.list <- c(ri.list, rit) # Store the retention index in the list
      print(paste("RI:", rit, sep = " "))  # Report RI
    }
    Lap = Lap + 1 # adjust the lap count
  }
  
  #remove standard peaks
  Comps2 <- Comps2 %>% filter(!Compound %in% (Comps %>% 
                                                filter(Class == "STD") %>% 
                                                pull(Compound)))
  
  # Add RI columns with all the calculated RI to the Comps data frame
  # replacing the mean_RT column
  Comps2 <- cbind.data.frame(Comps2 %>% select(-mean_RT), "RI" = ri.list)
  
  cat('\n')
  print("The RI calculation of all compounds has finished!") 
  cat('\n')
  
  # Return the new Comps data frame
  Comps2 %>% as_tibble()
}

# Import data files ----
## Load aligned GC-MS tables generated by script 02
load(here("AMelliMelli-CHC-data", "processed", "aligned_gcms-data.Rdata"))

## Load the compound information, result from the identification process
path_comps_id <- list.files(path = here("AMelliMelli-CHC-data", "tmp")
                            , pattern = "compounds-id.csv"
                            , full.names = T)
comps_id <- lapply(path_comps_id, read_csv)
summary(comps_id)

names(comps_id) <- str_split(path_comps_id
                                 , "/"
                                 , simplify = T) %>% 
  str_subset(".csv") %>% 
  str_remove("_compounds-id.csv") 
summary(comps_id)
comps_id

## Load group membership information
grouping_info <- read_csv(here("AMelliMelli-CHC-data"
                               , "raw"
                               , "AMelliMelli_samples-list.csv"))
grouping_info

# Shape the data frames ----
## Standards ####
### Extract standards information from comps_id
std_info <- cbind(comps_id$std %>% 
                    select(Compound, mean_RT)
                  , data.frame(row.names = rownames(comps_id$std)
                               , t(as.data.frame(strsplit(comps_id$std$Compound
                                                          , "_"))))) %>% 
  as_tibble()

### Define columns names
colnames(std_info) <- c("Compound", "mean_RT", "Chain.length"
                        , "Class", "Mod.position")
std_info

### Correct entries format
#### The peaks/compounds within the standard runs, must be differentiated from
#### the alkanes in the samples
#### Class
std_info['Class'][std_info['Class'] == "ane"] <- "STD"
std_info

#### Compound names
# Store the compound names in a vector, where they will be altered into their 
# final format
std_names <- std_info$Compound
# Set the name of the compound inside the standard samples
std_names[!is.na(std_names)] <- paste(std_info %>%
                                        filter(!is.na(Compound)) %>% 
                                        pull(Class)
                                      , std_info %>%
                                        filter(!is.na(Compound)) %>% 
                                        pull(Chain.length)
                                      , sep = "-")
# Change the compound names in the std_info data frame to the correct format names
std_info$Compound <- std_names
std_info

#### Chain length
# It needs to be defined as an integer
std_info$Chain.length <- std_info$Chain.length %>%
  str_remove("C") %>%
  as.integer()
std_info

## Compounds information ####
### IW ####
### Extract compounds information from comps_id
IW_comps_info <- cbind(comps_id$IW
                  , data.frame(row.names = rownames(comps_id$IW)
                               , t(as.data.frame(strsplit(comps_id$IW$Compound
                                                          , "_"))))) %>% 
  as_tibble()

### Define columns names
colnames(IW_comps_info) <- c(colnames(comps_id$IW)
                          , "Chain.length", "Class", "Mod.position")
IW_comps_info

### Correct entries format
#### Compound names
# Store the compound names in a vector, where they will be altered into their final format
comps_names <- IW_comps_info$Compound

# Set the first part of the identified CHC names (full name in the case of the alkanes)
comps_names[!is.na(comps_names)] <- paste0(IW_comps_info %>% 
                                             filter(!is.na(Compound)) %>% 
                                             pull(Mod.position)
                                           , IW_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Class)
                                           , IW_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Chain.length)) %>% 
  str_replace_all(paste(c("ane", "ene", "diene"), collapse = '|'), "NA") %>% 
  str_remove_all("NA")

# Finish formatting the names of the identified CHC
comps_names[!is.na(comps_names)] <- paste(comps_names[!is.na(comps_names)]
                                          , IW_comps_info %>% 
                                            filter(!is.na(Compound)) %>% 
                                            pull(Class) %>% 
                                            str_replace(paste(c("ane", "Me"
                                                                , "Dime", "Trime")
                                                              , collapse = '|')
                                                        , "NA")
                                          , sep = ":") %>% 
  str_remove(":NA")

# Change the compound names in the comps_info data frame to the correct format names
IW_comps_info$Compound <- comps_names
IW_comps_info
IW_comps_info %>% filter(!is.na(Compound))

#### Chain length
IW_comps_info$Chain.length <- IW_comps_info$Chain.length %>%
  str_remove("C") %>%
  as.integer()
IW_comps_info
IW_comps_info %>% filter(!is.na(Compound))

#### Class
IW_comps_info['Class'][IW_comps_info['Class'] == "ane"] <- "Alkane"
IW_comps_info['Class'][IW_comps_info['Class'] == "ene"] <- "Alkene"
IW_comps_info['Class'][IW_comps_info['Class'] == "diene"] <- "Alkadiene"
IW_comps_info['Class'][IW_comps_info['Class'] == "Me"] <- "Methyl"
IW_comps_info['Class'][IW_comps_info['Class'] == "Dime"] <- "Dimethyl"
IW_comps_info['Class'][IW_comps_info['Class'] == "Trime"] <- "Trimethyl"
IW_comps_info['Class'][IW_comps_info['Class'] == "Tetrame"] <- "Tetramethyl"
IW_comps_info
IW_comps_info %>% filter(!is.na(Compound))

#### Add mean_RT 
IW_comps_info <- my_merge(IW_comps_info
                       , IW_RT %>% 
                         select(Peak, mean_RT)
                       , merge_by = "Peak") %>% 
  arrange(mean_RT) %>% 
  as_tibble()
IW_comps_info
IW_comps_info %>% filter(!is.na(Compound))

### OW ####
### Extract compounds information from comps_id
OW_comps_info <- cbind(comps_id$OW
                       , data.frame(row.names = rownames(comps_id$OW)
                                    , t(as.data.frame(strsplit(comps_id$OW$Compound
                                                               , "_"))))) %>% 
  as_tibble()

### Define columns names
colnames(OW_comps_info) <- c(colnames(comps_id$OW)
                             , "Chain.length", "Class", "Mod.position")
OW_comps_info

### Correct entries format
#### Compound names
# Store the compound names in a vector, where they will be altered into their final format
comps_names <- OW_comps_info$Compound

# Set the first part of the identified CHC names (full name in the case of the alkanes)
comps_names[!is.na(comps_names)] <- paste0(OW_comps_info %>% 
                                             filter(!is.na(Compound)) %>% 
                                             pull(Mod.position)
                                           , OW_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Class)
                                           , OW_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Chain.length)) %>% 
  str_replace_all(paste(c("ane", "ene", "diene"), collapse = '|'), "NA") %>% 
  str_remove_all("NA")

# Finish formatting the names of the identified CHC
comps_names[!is.na(comps_names)] <- paste(comps_names[!is.na(comps_names)]
                                          , OW_comps_info %>% 
                                            filter(!is.na(Compound)) %>% 
                                            pull(Class) %>% 
                                            str_replace(paste(c("ane", "Me")
                                                              , collapse = '|')
                                                        , "NA")
                                          , sep = ":") %>% 
  str_remove(":NA")

# Change the compound names in the comps_info data frame to the correct format names
OW_comps_info$Compound <- comps_names
OW_comps_info
OW_comps_info %>% filter(!is.na(Compound))

#### Chain length
OW_comps_info$Chain.length <- OW_comps_info$Chain.length %>%
  str_remove("C") %>%
  as.integer()
OW_comps_info
OW_comps_info %>% filter(!is.na(Compound))

#### Class
OW_comps_info['Class'][OW_comps_info['Class'] == "ane"] <- "Alkane"
OW_comps_info['Class'][OW_comps_info['Class'] == "ene"] <- "Alkene"
OW_comps_info['Class'][OW_comps_info['Class'] == "diene"] <- "Alkadiene"
OW_comps_info['Class'][OW_comps_info['Class'] == "Me"] <- "Methyl"
OW_comps_info['Class'][OW_comps_info['Class'] == "Dime"] <- "Dimethyl"
OW_comps_info['Class'][OW_comps_info['Class'] == "Trime"] <- "Trimethyl"
OW_comps_info['Class'][OW_comps_info['Class'] == "Tetrame"] <- "Tetramethyl"
OW_comps_info
OW_comps_info %>% filter(!is.na(Compound))

#### Add mean_RT 
OW_comps_info <- my_merge(OW_comps_info
                          , OW_RT %>% 
                            select(Peak, mean_RT)
                          , merge_by = "Peak") %>% 
  arrange(mean_RT) %>% 
  as_tibble()
OW_comps_info
OW_comps_info %>% filter(!is.na(Compound))

## Samples ####
### IW ####
### Add CHC identification to the aligned sample data frames
IW_RT <- my_merge(IW_comps_info %>% 
                    select(Peak, Compound)
                  , IW_RT %>% 
                    select(!mean_RT)
                  , merge_by = "Peak") %>%
  as_tibble()
IW_RT

IW_area <- my_merge(IW_comps_info %>% 
                      select(Peak, Compound)
                    , IW_area %>% 
                      select(!mean_RT)
                    , merge_by = "Peak") %>%
  as_tibble()
IW_area

### OW ####
### Add CHC identification to the aligned sample data frames
OW_RT <- my_merge(OW_comps_info %>% 
                    select(Peak, Compound)
                  , OW_RT %>% 
                    select(!mean_RT)
                  , merge_by = "Peak") %>%
  as_tibble()
OW_RT

OW_area <- my_merge(OW_comps_info %>% 
                      select(Peak, Compound)
                    , OW_area %>% 
                      select(!mean_RT)
                    , merge_by = "Peak") %>%
  as_tibble()
OW_area

### grouping_info ####
### Define grouping variables as factors
grouping_info[, colnames(grouping_info)] <-
  lapply(grouping_info[, colnames(grouping_info)], as.factor)
grouping_info

# Clean the data ----
## Delete trace compounds ####
### Extract the abundance data into a data frame
IW_daten <- IW_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

OW_daten <- OW_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

### Use trace_comps function to delete compounds, from each sample, 
### with an abundance below a defined threshold
IW_daten <- trace_comps(IW_daten, threshold = 0.01)
OW_daten <- trace_comps(OW_daten, threshold = 0.01)

### Adjust the area data frame
IW_area <- cbind(IW_area %>% 
                   select(Peak:Compound)
                 , IW_daten) %>% 
  as_tibble()
IW_area

OW_area <- cbind(OW_area %>% 
                   select(Peak:Compound)
                 , OW_daten) %>% 
  as_tibble()
OW_area

### Adjust the RT data frame accordingly
IW_RT_daten <- IW_RT %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

OW_RT_daten <- OW_RT %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

IW_RT_daten[is.na(IW_daten)] <- NA
IW_RT <- cbind(IW_RT %>% 
                 select(Peak:Compound)
               , IW_RT_daten) %>% 
  as_tibble()
IW_RT

OW_RT_daten[is.na(OW_daten)] <- NA
OW_RT <- cbind(OW_RT %>% 
                 select(Peak:Compound)
               , OW_RT_daten) %>% 
  as_tibble()
OW_RT

## Delete non-CHC compounds ####
### Samples 
IW_area <- IW_area %>% filter(!is.na(Compound))
IW_area
IW_RT <- IW_RT %>% filter(!is.na(Compound))
IW_RT

OW_area <- OW_area %>% filter(!is.na(Compound))
OW_area
OW_RT <- OW_RT %>% filter(!is.na(Compound))
OW_RT

### Compounds information
IW_comps_info <- IW_comps_info %>% filter(!is.na(Compound))
IW_comps_info

OW_comps_info <- OW_comps_info %>% filter(!is.na(Compound))
OW_comps_info

### Standards
std_info <- std_info %>% filter(!is.na(Compound))
std_info

## Delete rare compounds - OPTIONAL ####
### IW ####
IW_area <- below_50(IW_area
                    , grouping.info = grouping_info %>% 
                      filter(Task == "In-hive workers"))
IW_area

IW_RT <- IW_RT %>% 
  filter(Peak %in% IW_area$Peak)
IW_RT

IW_comps_info <- IW_comps_info %>% 
  filter(Peak %in% IW_area$Peak)
IW_comps_info

### OW ####
OW_area <- below_50(OW_area
                    , grouping.info = grouping_info %>% 
                      filter(Task == "Out-hive workers"))
OW_area

OW_RT <- OW_RT %>% 
  filter(Peak %in% OW_area$Peak)
OW_RT

OW_comps_info <- OW_comps_info %>% 
  filter(Peak %in% OW_area$Peak)
OW_comps_info

# Kováts Retention index ----
## Recalculate mean_RT ####
# It may have changed due to the fusion and deletion of peaks
# across the data set  
### IW ####
identical(IW_RT %>%
            select(!Peak:Compound) %>%
            rowMeans(na.rm = T)
          , IW_comps_info$mean_RT)

IW_comps_info$mean_RT <- IW_RT %>%
  select(!Peak:Compound) %>% 
  rowMeans(na.rm = T)

identical(IW_RT %>%
            select(!Peak:Compound) %>% 
            rowMeans(na.rm = T)
          , IW_comps_info$mean_RT)
IW_comps_info$mean_RT
IW_comps_info

### OW ####
identical(OW_RT %>%
            select(!Peak:Compound) %>%
            rowMeans(na.rm = T)
          , OW_comps_info$mean_RT)

OW_comps_info$mean_RT <- OW_RT %>%
  select(!Peak:Compound) %>% 
  rowMeans(na.rm = T)

identical(OW_RT %>%
            select(!Peak:Compound) %>% 
            rowMeans(na.rm = T)
          , OW_comps_info$mean_RT)
OW_comps_info$mean_RT
OW_comps_info

## Delete unnecessary standards ####
IW_std_info <- std_info %>% 
  filter(!Chain.length %in% (IW_comps_info %>% 
           filter(Class == "Alkane") %>% 
           pull(Chain.length)))
IW_std_info

OW_std_info <- std_info %>% 
  filter(!Chain.length %in% (OW_comps_info %>% 
                               filter(Class == "Alkane") %>% 
                               pull(Chain.length)))
OW_std_info

## Calculate RI ####
# Merge comps_info and std_info, so the mean_Rt of every alkane 
# is in the same data frame
IW_comps_info <- merge(IW_std_info
                       , IW_comps_info
                       , all = T
                       , sort = F) %>%
  arrange(mean_RT) %>% 
  select(Peak
         , everything()) %>% 
  as_tibble()
IW_comps_info

OW_comps_info <- merge(OW_std_info
                       , OW_comps_info
                       , all = T
                       , sort = F) %>%
  arrange(mean_RT) %>% 
  select(Peak
         , everything()) %>% 
  as_tibble()
OW_comps_info

# Use the k.ri() function to calculate the retention index
IW_comps_info <- k.ri(IW_comps_info)
OW_comps_info <- k.ri(OW_comps_info)

# Relative abundance  (%) ----
## Extract the abundance data into a data frame
IW_daten <- IW_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

OW_daten <- OW_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

## Replace NAs with 0s
IW_daten[is.na(IW_daten)] <- 0
OW_daten[is.na(OW_daten)] <- 0

### Calculate the relative abundance (%) of each compound per sample 
IW_daten <- IW_daten %>% t / rowSums(IW_daten %>% t) * 100
IW_daten <- IW_daten %>% 
  t %>% 
  as.data.frame() %>% 
  as_tibble()

OW_daten <- OW_daten %>% t / rowSums(OW_daten %>% t) * 100
OW_daten <- OW_daten %>% 
  t %>% 
  as.data.frame() %>% 
  as_tibble()

### Verify that the sum of all relative abundances per sample is exactly 100
IW_daten %>% colSums()
OW_daten %>% colSums()

# Export final data frames ----
## Sort individuals' data columns according to the order of the individuals
## in the grouping_info data frames
IW_daten <- IW_daten %>% select(all_of(grouping_info %>% 
                                         filter(Task == "In-hive workers") %>% 
                                         pull(Individual)))

OW_daten <- OW_daten %>% select(all_of(grouping_info %>% 
                                         filter(Task == "Out-hive workers") %>% 
                                         pull(Individual)))

## Adjust the area data frame
IW_area <- cbind(IW_area %>% 
                   select(Peak:Compound)
                 , IW_daten) %>% 
  as_tibble()

OW_area <- cbind(OW_area %>% 
                   select(Peak:Compound)
                 , OW_daten) %>% 
  as_tibble()

## merge the area data frame and comps_info data frame to adjust peak numbering
IW_table <- my_merge(IW_comps_info
                  , IW_area
                  , merge_by = c("Peak", "Compound")) %>% 
  as.data.frame()
IW_table
str(IW_table)

OW_table <- my_merge(OW_comps_info
                  , OW_area
                  , merge_by = c("Peak", "Compound")) %>% 
  as.data.frame()
OW_table
str(OW_table)

## Adjust peak numbering
IW_table$Peak <- paste0("P", 1:length(IW_table$Peak))
head(IW_table)
str(IW_table)

OW_table$Peak <- paste0("P", 1:length(OW_table$Peak))
head(OW_table)
str(OW_table)

## Define peak numbers as row names
rownames(IW_table) <- IW_table$Peak
IW_table <- IW_table %>% select(-Peak)
head(IW_table)
str(IW_table)

rownames(OW_table) <- OW_table$Peak
OW_table <- OW_table %>% select(-Peak)
head(OW_table)
str(OW_table)

# ## re-split composition data frame and comps_info
# IW_daten <- IW_table %>% 
#   select(!Compound:RI)
# head(IW_daten)
# str(IW_daten)
# 
# IW_comps_info <- IW_table %>% select(Compound:RI)
# head(IW_comps_info)
# str(IW_comps_info)
# 
# OW_daten <- OW_table %>% 
#   select(!Compound:RI)
# head(OW_daten)
# str(OW_daten)
# 
# OW_comps_info <- OW_table %>% select(Compound:RI)
# head(OW_comps_info)
# str(OW_comps_info)

## Export the data frames
save(list = c("IW_table"
              , "OW_table"
              , "grouping_info")
     , file = here("AMelliMelli-CHC-data"
                   , "processed"
                   , "group_tables.Rdata"))
print("The data frames for analysis were exported")

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "SInf_Script03.txt"))
print("The sessionInfo report was exported. The script 03 finished running")

# ## Detach/unload packages
# lapply(NPacks, unloadNamespace)
# 
# ## Clear environment
# rm(list=ls())
