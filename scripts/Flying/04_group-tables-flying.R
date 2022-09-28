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

## Modified read CSV function for comps_id
read_csv_comps_id <- function(file.path) {
  read_csv(file = file.path, col_select = c("Peak", "Compound", "mean_RT"))
  
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
    unite('group', !c("Bee_number", "Task", "Hive"), remove = F)
  
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
                  pull(Bee_number) %>% 
                  as.character() %>% 
                  paste0(collapse = " ")
                , sep = " "))
    
    df_g <- df_area %>% 
      select(all_of(grouping.info %>% 
                      filter(group == g) %>% 
                      pull(Bee_number)))
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
  levels(Comps2$Class) <- c(levels(Comps2$Class), "Alkane")
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
  print(Comps2)
  cat('\n')
  print("The RI calculation of all compounds has finished!") 
  cat('\n')
  
  # Return the new Comps data frame
  Comps2 %>% as_tibble()
}

# Import data files ----
## Load aligned GC-MS tables generated by script 03
load(here("data", "processed", "Flying", "aligned_gcms-data.Rdata"))

## Load the compound information, result from the identification process
path_comps_id <- list.files(path = here("data", "raw", "Flying", "tmp")
                            , pattern = "ID_table.csv"
                            , full.names = T)
comps_id <- lapply(path_comps_id, read_csv_comps_id)
summary(comps_id)

names(comps_id) <- str_split(path_comps_id
                                 , "/"
                                 , simplify = T) %>% 
  str_subset(".csv") %>% 
  str_remove("_RT_ID_table.csv") 
summary(comps_id)
comps_id

## Load group membership information
grouping_info <- read_csv(here("data"
                               , "raw"
                               , "Flying"
                               , "Samples-Flying.csv"))
grouping_info

# Shape the data frames ----
## Standards ####
### Extract standards information from comps_id
std_info <- cbind(comps_id$STD %>% 
                    select(Compound, mean_RT)
                  , data.frame(row.names = rownames(comps_id$STD)
                               , t(as.data.frame(strsplit(comps_id$STD$Compound
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
levels(std_info$Class) <- c(levels(std_info$Class), "STD")
std_info$Class[std_info$Class == "ane"] <- "STD"
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
### Ca_false ####
### Extract compounds information from comps_id
Ca_false_comps_info <- cbind(comps_id$Ca_false
                  , data.frame(row.names = rownames(comps_id$Ca_false)
                               , t(as.data.frame(strsplit(comps_id$Ca_false$Compound
                                                          , "_"))))) %>% 
  as_tibble()

### Define columns names
colnames(Ca_false_comps_info) <- c(colnames(comps_id$Ca_false)
                          , "Chain.length", "Class", "Mod.position")
Ca_false_comps_info

### Correct entries format
#### Compound names
# Store the compound names in a vector, where they will be altered into their final format
comps_names <- Ca_false_comps_info$Compound

# Set the first part of the identified CHC names (full name in the case of the alkanes)
comps_names[!is.na(comps_names)] <- paste0(Ca_false_comps_info %>% 
                                             filter(!is.na(Compound)) %>% 
                                             pull(Mod.position)
                                           , Ca_false_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Class)
                                           , Ca_false_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Chain.length)) %>% 
  str_replace_all(paste(c("ane", "ene", "diene"), collapse = '|'), "NA") %>% 
  str_remove_all("NA")

# Finish formatting the names of the identified CHC
comps_names[!is.na(comps_names)] <- paste(comps_names[!is.na(comps_names)]
                                          , Ca_false_comps_info %>% 
                                            filter(!is.na(Compound)) %>% 
                                            pull(Class) %>% 
                                            str_replace(paste(c("ane", "Me"
                                                                , "Dime", "Trime")
                                                              , collapse = '|')
                                                        , "NA")
                                          , sep = ":") %>% 
  str_remove(":NA")

# Change the compound names in the comps_info data frame to the correct format names
Ca_false_comps_info$Compound <- comps_names
Ca_false_comps_info
Ca_false_comps_info %>% filter(!is.na(Compound))

#### Chain length
Ca_false_comps_info$Chain.length <- Ca_false_comps_info$Chain.length %>%
  str_remove("C") %>%
  as.integer()
Ca_false_comps_info
Ca_false_comps_info %>% filter(!is.na(Compound))

#### Class
levels(Ca_false_comps_info$Class) <- c(levels(Ca_false_comps_info$Class)
                                       , "Alkane"
                                       , "Alkene"
                                       , "Alkadiene"
                                       , "Methyl"
                                       , "Dimethyl"
                                       , "Trimethyl"
                                       , "Tetramethyl")

Ca_false_comps_info['Class'][Ca_false_comps_info['Class'] == "ane"] <- "Alkane"
Ca_false_comps_info['Class'][Ca_false_comps_info['Class'] == "ene"] <- "Alkene"
Ca_false_comps_info['Class'][Ca_false_comps_info['Class'] == "diene"] <- "Alkadiene"
Ca_false_comps_info['Class'][Ca_false_comps_info['Class'] == "Me"] <- "Methyl"
Ca_false_comps_info['Class'][Ca_false_comps_info['Class'] == "Dime"] <- "Dimethyl"
Ca_false_comps_info['Class'][Ca_false_comps_info['Class'] == "Trime"] <- "Trimethyl"
Ca_false_comps_info['Class'][Ca_false_comps_info['Class'] == "Tetrame"] <- "Tetramethyl"
Ca_false_comps_info
Ca_false_comps_info %>% filter(!is.na(Compound))

#### Add mean_RT 
# Ca_false_RT <- mg_list_RT[["Ca_FALSE"]]
# 
# Ca_false_comps_info <- my_merge(Ca_false_comps_info
#                        , Ca_false_RT %>% 
#                          select(Peak, mean_RT)
#                        , merge_by = "Peak") %>% 
#   arrange(mean_RT) %>% 
#   as_tibble()
# Ca_false_comps_info
# Ca_false_comps_info %>% filter(!is.na(Compound))

### Ca_true ####
### Extract compounds information from comps_id
Ca_true_comps_info <- cbind(comps_id$Ca_true
                       , data.frame(row.names = rownames(comps_id$Ca_true)
                                    , t(as.data.frame(strsplit(comps_id$Ca_true$Compound
                                                               , "_"))))) %>% 
  as_tibble()

### Define columns names
colnames(Ca_true_comps_info) <- c(colnames(comps_id$Ca_true)
                             , "Chain.length", "Class", "Mod.position")
Ca_true_comps_info

### Correct entries format
#### Compound names
# Store the compound names in a vector, where they will be altered into their final format
comps_names <- Ca_true_comps_info$Compound

# Set the first part of the identified CHC names (full name in the case of the alkanes)
comps_names[!is.na(comps_names)] <- paste0(Ca_true_comps_info %>% 
                                             filter(!is.na(Compound)) %>% 
                                             pull(Mod.position)
                                           , Ca_true_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Class)
                                           , Ca_true_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Chain.length)) %>% 
  str_replace_all(paste(c("ane", "ene", "diene"), collapse = '|'), "NA") %>% 
  str_remove_all("NA")

# Finish formatting the names of the identified CHC
comps_names[!is.na(comps_names)] <- paste(comps_names[!is.na(comps_names)]
                                          , Ca_true_comps_info %>% 
                                            filter(!is.na(Compound)) %>% 
                                            pull(Class) %>% 
                                            str_replace(paste(c("ane", "Me")
                                                              , collapse = '|')
                                                        , "NA")
                                          , sep = ":") %>% 
  str_remove(":NA")

# Change the compound names in the comps_info data frame to the correct format names
Ca_true_comps_info$Compound <- comps_names
Ca_true_comps_info
Ca_true_comps_info %>% filter(!is.na(Compound))

#### Chain length
Ca_true_comps_info$Chain.length <- Ca_true_comps_info$Chain.length %>%
  str_remove("C") %>%
  as.integer()
Ca_true_comps_info
Ca_true_comps_info %>% filter(!is.na(Compound))

#### Class
levels(Ca_true_comps_info$Class) <- c(levels(Ca_true_comps_info$Class)
                                       , "Alkane"
                                       , "Alkene"
                                       , "Alkadiene"
                                       , "Methyl"
                                       , "Dimethyl"
                                       , "Trimethyl"
                                       , "Tetramethyl")

Ca_true_comps_info['Class'][Ca_true_comps_info['Class'] == "ane"] <- "Alkane"
Ca_true_comps_info['Class'][Ca_true_comps_info['Class'] == "ene"] <- "Alkene"
Ca_true_comps_info['Class'][Ca_true_comps_info['Class'] == "diene"] <- "Alkadiene"
Ca_true_comps_info['Class'][Ca_true_comps_info['Class'] == "Me"] <- "Methyl"
Ca_true_comps_info['Class'][Ca_true_comps_info['Class'] == "Dime"] <- "Dimethyl"
Ca_true_comps_info['Class'][Ca_true_comps_info['Class'] == "Trime"] <- "Trimethyl"
Ca_true_comps_info['Class'][Ca_true_comps_info['Class'] == "Tetrame"] <- "Tetramethyl"
Ca_true_comps_info
Ca_true_comps_info %>% filter(!is.na(Compound))

# #### Add mean_RT 
# Ca_true_comps_info <- my_merge(Ca_true_comps_info
#                           , Ca_true_RT %>% 
#                             select(Peak, mean_RT)
#                           , merge_by = "Peak") %>% 
#   arrange(mean_RT) %>% 
#   as_tibble()
# Ca_true_comps_info
# Ca_true_comps_info %>% filter(!is.na(Compound))

### Ib_false ####
### Extract compounds information from comps_id
Ib_false_comps_info <- cbind(comps_id$Ib_false
                            , data.frame(row.names = rownames(comps_id$Ib_false)
                                         , t(as.data.frame(strsplit(comps_id$Ib_false$Compound
                                                                    , "_"))))) %>% 
  as_tibble()

### Define columns names
colnames(Ib_false_comps_info) <- c(colnames(comps_id$Ib_false)
                                  , "Chain.length", "Class", "Mod.position")
Ib_false_comps_info

### Correct entries format
#### Compound names
# Store the compound names in a vector, where they will be altered into their final format
comps_names <- Ib_false_comps_info$Compound

# Set the first part of the identified CHC names (full name in the case of the alkanes)
comps_names[!is.na(comps_names)] <- paste0(Ib_false_comps_info %>% 
                                             filter(!is.na(Compound)) %>% 
                                             pull(Mod.position)
                                           , Ib_false_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Class)
                                           , Ib_false_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Chain.length)) %>% 
  str_replace_all(paste(c("ane", "ene", "diene"), collapse = '|'), "NA") %>% 
  str_remove_all("NA")

# Finish formatting the names of the identified CHC
comps_names[!is.na(comps_names)] <- paste(comps_names[!is.na(comps_names)]
                                          , Ib_false_comps_info %>% 
                                            filter(!is.na(Compound)) %>% 
                                            pull(Class) %>% 
                                            str_replace(paste(c("ane", "Me")
                                                              , collapse = '|')
                                                        , "NA")
                                          , sep = ":") %>% 
  str_remove(":NA")

# Change the compound names in the comps_info data frame to the correct format names
Ib_false_comps_info$Compound <- comps_names
Ib_false_comps_info
Ib_false_comps_info %>% filter(!is.na(Compound))

#### Chain length
Ib_false_comps_info$Chain.length <- Ib_false_comps_info$Chain.length %>%
  str_remove("C") %>%
  as.integer()
Ib_false_comps_info
Ib_false_comps_info %>% filter(!is.na(Compound))

#### Class
levels(Ib_false_comps_info$Class) <- c(levels(Ib_false_comps_info$Class)
                                       , "Alkane"
                                       , "Alkene"
                                       , "Alkadiene"
                                       , "Methyl"
                                       , "Dimethyl"
                                       , "Trimethyl"
                                       , "Tetramethyl")

Ib_false_comps_info['Class'][Ib_false_comps_info['Class'] == "ane"] <- "Alkane"
Ib_false_comps_info['Class'][Ib_false_comps_info['Class'] == "ene"] <- "Alkene"
Ib_false_comps_info['Class'][Ib_false_comps_info['Class'] == "diene"] <- "Alkadiene"
Ib_false_comps_info['Class'][Ib_false_comps_info['Class'] == "Me"] <- "Methyl"
Ib_false_comps_info['Class'][Ib_false_comps_info['Class'] == "Dime"] <- "Dimethyl"
Ib_false_comps_info['Class'][Ib_false_comps_info['Class'] == "Trime"] <- "Trimethyl"
Ib_false_comps_info['Class'][Ib_false_comps_info['Class'] == "Tetrame"] <- "Tetramethyl"
Ib_false_comps_info
Ib_false_comps_info %>% filter(!is.na(Compound))

#### Add mean_RT 
# Ib_false_comps_info <- my_merge(Ib_false_comps_info
#                                , Ib_false_RT %>% 
#                                  select(Peak, mean_RT)
#                                , merge_by = "Peak") %>% 
#   arrange(mean_RT) %>% 
#   as_tibble()
# Ib_false_comps_info
# Ib_false_comps_info %>% filter(!is.na(Compound))

### Ib_true ####
### Extract compounds information from comps_id
Ib_true_comps_info <- cbind(comps_id$Ib_true
                            , data.frame(row.names = rownames(comps_id$Ib_true)
                                         , t(as.data.frame(strsplit(comps_id$Ib_true$Compound
                                                                    , "_"))))) %>% 
  as_tibble()

### Define columns names
colnames(Ib_true_comps_info) <- c(colnames(comps_id$Ib_true)
                                  , "Chain.length", "Class", "Mod.position")
Ib_true_comps_info

### Correct entries format
#### Compound names
# Store the compound names in a vector, where they will be altered into their final format
comps_names <- Ib_true_comps_info$Compound

# Set the first part of the identified CHC names (full name in the case of the alkanes)
comps_names[!is.na(comps_names)] <- paste0(Ib_true_comps_info %>% 
                                             filter(!is.na(Compound)) %>% 
                                             pull(Mod.position)
                                           , Ib_true_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Class)
                                           , Ib_true_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Chain.length)) %>% 
  str_replace_all(paste(c("ane", "ene", "diene"), collapse = '|'), "NA") %>% 
  str_remove_all("NA")

# Finish formatting the names of the identified CHC
comps_names[!is.na(comps_names)] <- paste(comps_names[!is.na(comps_names)]
                                          , Ib_true_comps_info %>% 
                                            filter(!is.na(Compound)) %>% 
                                            pull(Class) %>% 
                                            str_replace(paste(c("ane", "Me")
                                                              , collapse = '|')
                                                        , "NA")
                                          , sep = ":") %>% 
  str_remove(":NA")

# Change the compound names in the comps_info data frame to the correct format names
Ib_true_comps_info$Compound <- comps_names
Ib_true_comps_info
Ib_true_comps_info %>% filter(!is.na(Compound))

#### Chain length
Ib_true_comps_info$Chain.length <- Ib_true_comps_info$Chain.length %>%
  str_remove("C") %>%
  as.integer()
Ib_true_comps_info
Ib_true_comps_info %>% filter(!is.na(Compound))

#### Class
levels(Ib_true_comps_info$Class) <- c(levels(Ib_true_comps_info$Class)
                                       , "Alkane"
                                       , "Alkene"
                                       , "Alkadiene"
                                       , "Methyl"
                                       , "Dimethyl"
                                       , "Trimethyl"
                                       , "Tetramethyl")


Ib_true_comps_info['Class'][Ib_true_comps_info['Class'] == "ane"] <- "Alkane"
Ib_true_comps_info['Class'][Ib_true_comps_info['Class'] == "ene"] <- "Alkene"
Ib_true_comps_info['Class'][Ib_true_comps_info['Class'] == "diene"] <- "Alkadiene"
Ib_true_comps_info['Class'][Ib_true_comps_info['Class'] == "Me"] <- "Methyl"
Ib_true_comps_info['Class'][Ib_true_comps_info['Class'] == "Dime"] <- "Dimethyl"
Ib_true_comps_info['Class'][Ib_true_comps_info['Class'] == "Trime"] <- "Trimethyl"
Ib_true_comps_info['Class'][Ib_true_comps_info['Class'] == "Tetrame"] <- "Tetramethyl"
Ib_true_comps_info
Ib_true_comps_info %>% filter(!is.na(Compound))

# #### Add mean_RT 
# Ib_true_comps_info <- my_merge(Ib_true_comps_info
#                                , Ib_true_RT %>% 
#                                  select(Peak, mean_RT)
#                                , merge_by = "Peak") %>% 
#   arrange(mean_RT) %>% 
#   as_tibble()
# Ib_true_comps_info
# Ib_true_comps_info %>% filter(!is.na(Compound))

## Samples ####
### Ca_false ####
### Add CHC identification to the aligned sample data frames
Ca_false_RT <- my_merge(Ca_false_comps_info %>% 
                    select(Peak, Compound)
                  , mg_list_RT[["Ca_FALSE"]] %>% 
                    select(!mean_RT)
                  , merge_by = "Peak") %>%
  as_tibble()
Ca_false_RT

Ca_false_area <- my_merge(Ca_false_comps_info %>% 
                      select(Peak, Compound)
                    , mg_list_area[["Ca_FALSE"]] %>% 
                      select(!mean_RT)
                    , merge_by = "Peak") %>%
  as_tibble()
Ca_false_area

### Ca_true ####
### Add CHC identification to the aligned sample data frames
Ca_true_RT <- my_merge(Ca_true_comps_info %>% 
                    select(Peak, Compound)
                  , mg_list_RT[["Ca_TRUE"]] %>% 
                    select(!mean_RT)
                  , merge_by = "Peak") %>%
  as_tibble()
Ca_true_RT

Ca_true_area <- my_merge(Ca_true_comps_info %>% 
                      select(Peak, Compound)
                    , mg_list_area[["Ca_TRUE"]] %>% 
                      select(!mean_RT)
                    , merge_by = "Peak") %>%
  as_tibble()
Ca_true_area

### Ib_false ####
### Add CHC identification to the aligned sample data frames
Ib_false_RT <- my_merge(Ib_false_comps_info %>% 
                         select(Peak, Compound)
                       , mg_list_RT[["Ib_FALSE"]] %>% 
                         select(!mean_RT)
                       , merge_by = "Peak") %>%
  as_tibble()
Ib_false_RT

Ib_false_area <- my_merge(Ib_false_comps_info %>% 
                           select(Peak, Compound)
                         , mg_list_area[["Ib_FALSE"]] %>% 
                           select(!mean_RT)
                         , merge_by = "Peak") %>%
  as_tibble()
Ib_false_area

### Ib_true ####
### Add CHC identification to the aligned sample data frames
Ib_true_RT <- my_merge(Ib_true_comps_info %>% 
                         select(Peak, Compound)
                       , mg_list_RT[["Ib_TRUE"]] %>% 
                         select(!mean_RT)
                       , merge_by = "Peak") %>%
  as_tibble()
Ib_true_RT

Ib_true_area <- my_merge(Ib_true_comps_info %>% 
                           select(Peak, Compound)
                         , mg_list_area[["Ib_TRUE"]] %>% 
                           select(!mean_RT)
                         , merge_by = "Peak") %>%
  as_tibble()
Ib_true_area

### grouping_info ####
### Define grouping variables as factors
grouping_info[, colnames(grouping_info)] <-
  lapply(grouping_info[, colnames(grouping_info)], as.factor)
grouping_info

# Clean the data ----
## Delete trace compounds ####
### Extract the abundance data into a data frame
Ca_false_daten <- Ca_false_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ca_true_daten <- Ca_true_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ib_false_daten <- Ib_false_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ib_true_daten <- Ib_true_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

### Use trace_comps function to delete compounds, from each sample, 
### with an abundance below a defined threshold
Ca_false_daten <- trace_comps(Ca_false_daten, threshold = 0.01)
Ca_true_daten <- trace_comps(Ca_true_daten, threshold = 0.01)
Ib_false_daten <- trace_comps(Ib_false_daten, threshold = 0.01)
Ib_true_daten <- trace_comps(Ib_true_daten, threshold = 0.01)

### Adjust the area data frame
Ca_false_area <- cbind(Ca_false_area %>% 
                   select(Peak:Compound)
                 , Ca_false_daten) %>% 
  as_tibble()
Ca_false_area

Ca_true_area <- cbind(Ca_true_area %>% 
                   select(Peak:Compound)
                 , Ca_true_daten) %>% 
  as_tibble()
Ca_true_area

Ib_false_area <- cbind(Ib_false_area %>% 
                        select(Peak:Compound)
                      , Ib_false_daten) %>% 
  as_tibble()
Ib_false_area

Ib_true_area <- cbind(Ib_true_area %>% 
                        select(Peak:Compound)
                      , Ib_true_daten) %>% 
  as_tibble()
Ib_true_area

### Adjust the RT data frame accordingly
Ca_false_RT_daten <- Ca_false_RT %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ca_true_RT_daten <- Ca_true_RT %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ib_false_RT_daten <- Ib_false_RT %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ib_true_RT_daten <- Ib_true_RT %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ca_false_RT_daten[is.na(Ca_false_daten)] <- NA
Ca_false_RT <- cbind(Ca_false_RT %>% 
                 select(Peak:Compound)
               , Ca_false_RT_daten) %>% 
  as_tibble()
Ca_false_RT

Ca_true_RT_daten[is.na(Ca_true_daten)] <- NA
Ca_true_RT <- cbind(Ca_true_RT %>% 
                 select(Peak:Compound)
               , Ca_true_RT_daten) %>% 
  as_tibble()
Ca_true_RT

Ib_false_RT_daten[is.na(Ib_false_daten)] <- NA
Ib_false_RT <- cbind(Ib_false_RT %>% 
                      select(Peak:Compound)
                    , Ib_false_RT_daten) %>% 
  as_tibble()
Ib_false_RT

Ib_true_RT_daten[is.na(Ib_true_daten)] <- NA
Ib_true_RT <- cbind(Ib_true_RT %>% 
                      select(Peak:Compound)
                    , Ib_true_RT_daten) %>% 
  as_tibble()
Ib_true_RT

## Delete non-CHC compounds ####
### Samples 
Ca_false_area <- Ca_false_area %>% filter(!is.na(Compound))
Ca_false_area
Ca_false_RT <- Ca_false_RT %>% filter(!is.na(Compound))
Ca_false_RT

Ca_true_area <- Ca_true_area %>% filter(!is.na(Compound))
Ca_true_area
Ca_true_RT <- Ca_true_RT %>% filter(!is.na(Compound))
Ca_true_RT

Ib_false_area <- Ib_false_area %>% filter(!is.na(Compound))
Ib_false_area
Ib_false_RT <- Ib_false_RT %>% filter(!is.na(Compound))
Ib_false_RT

Ib_true_area <- Ib_true_area %>% filter(!is.na(Compound))
Ib_true_area
Ib_true_RT <- Ib_true_RT %>% filter(!is.na(Compound))
Ib_true_RT

### Compounds information
Ca_false_comps_info <- Ca_false_comps_info %>% filter(!is.na(Compound))
Ca_false_comps_info

Ca_true_comps_info <- Ca_true_comps_info %>% filter(!is.na(Compound))
Ca_true_comps_info

Ib_true_comps_info <- Ib_true_comps_info %>% filter(!is.na(Compound))
Ib_true_comps_info

Ib_false_comps_info <- Ib_false_comps_info %>% filter(!is.na(Compound))
Ib_false_comps_info

### Standards
std_info <- std_info %>% filter(!is.na(Compound))
std_info

## Delete rare compounds - OPTIONAL ####
### Ca_false ####
Ca_false_area <- below_50(Ca_false_area
                    , grouping.info = grouping_info %>% 
                      filter(Subspecies == "Ca" 
                             & Flying == "FALSE"))
Ca_false_area

Ca_false_RT <- Ca_false_RT %>% 
  filter(Peak %in% Ca_false_area$Peak)
Ca_false_RT

Ca_false_comps_info <- Ca_false_comps_info %>% 
  filter(Peak %in% Ca_false_area$Peak)
Ca_false_comps_info

### Ca_true ####
Ca_true_area <- below_50(Ca_true_area
                    , grouping.info = grouping_info %>% 
                      filter(Subspecies == "Ca" 
                             & Flying == "TRUE"))
Ca_true_area

Ca_true_RT <- Ca_true_RT %>% 
  filter(Peak %in% Ca_true_area$Peak)
Ca_true_RT

Ca_true_comps_info <- Ca_true_comps_info %>% 
  filter(Peak %in% Ca_true_area$Peak)
Ca_true_comps_info

### Ib_false ####
Ib_false_area <- below_50(Ib_false_area
                              , grouping.info = grouping_info %>% 
                                filter(Subspecies == "Ib" 
                                       & Flying == "FALSE"))
Ib_false_area

Ib_false_RT <- Ib_false_RT %>% 
  filter(Peak %in% Ib_false_area$Peak)
Ib_false_RT

Ib_false_comps_info <- Ib_false_comps_info %>% 
  filter(Peak %in% Ib_false_area$Peak)
Ib_false_comps_info

### Ib_true ####
Ib_true_area <- below_50(Ib_true_area
                              , grouping.info = grouping_info %>% 
                                filter(Subspecies == "Ib" 
                                       & Flying == "TRUE"))
Ib_true_area

Ib_true_RT <- Ib_true_RT %>% 
  filter(Peak %in% Ib_true_area$Peak)
Ib_true_RT

Ib_true_comps_info <- Ib_true_comps_info %>% 
  filter(Peak %in% Ib_true_area$Peak)
Ib_true_comps_info

# Kováts Retention index ----
## Recalculate mean_RT ####
# It may have changed due to the fusion and deletion of peaks
# across the data set  
### Ca_false ####
identical(Ca_false_RT %>%
            select(!Peak:Compound) %>%
            rowMeans(na.rm = T)
          , Ca_false_comps_info$mean_RT)

Ca_false_comps_info$mean_RT <- Ca_false_RT %>%
  select(!Peak:Compound) %>% 
  rowMeans(na.rm = T)

identical(Ca_false_RT %>%
            select(!Peak:Compound) %>% 
            rowMeans(na.rm = T)
          , Ca_false_comps_info$mean_RT)
Ca_false_comps_info$mean_RT
Ca_false_comps_info

### Ca_true ####
identical(Ca_true_RT %>%
            select(!Peak:Compound) %>%
            rowMeans(na.rm = T)
          , Ca_true_comps_info$mean_RT)

Ca_true_comps_info$mean_RT <- Ca_true_RT %>%
  select(!Peak:Compound) %>% 
  rowMeans(na.rm = T)

identical(Ca_true_RT %>%
            select(!Peak:Compound) %>% 
            rowMeans(na.rm = T)
          , Ca_true_comps_info$mean_RT)
Ca_true_comps_info$mean_RT
Ca_true_comps_info

### Ib_false ####
identical(Ib_false_RT %>%
            select(!Peak:Compound) %>%
            rowMeans(na.rm = T)
          , Ib_false_comps_info$mean_RT)

Ib_false_comps_info$mean_RT <- Ib_false_RT %>%
  select(!Peak:Compound) %>% 
  rowMeans(na.rm = T)

identical(Ib_false_RT %>%
            select(!Peak:Compound) %>% 
            rowMeans(na.rm = T)
          , Ib_false_comps_info$mean_RT)
Ib_false_comps_info$mean_RT
Ib_false_comps_info


### Ib_true ####
identical(Ib_true_RT %>%
            select(!Peak:Compound) %>%
            rowMeans(na.rm = T)
          , Ib_true_comps_info$mean_RT)

Ib_true_comps_info$mean_RT <- Ib_true_RT %>%
  select(!Peak:Compound) %>% 
  rowMeans(na.rm = T)

identical(Ib_true_RT %>%
            select(!Peak:Compound) %>% 
            rowMeans(na.rm = T)
          , Ib_true_comps_info$mean_RT)
Ib_true_comps_info$mean_RT
Ib_true_comps_info

## Delete unnecessary standards ####
Ca_false_std_info <- std_info %>% 
  filter(!Chain.length %in% (Ca_false_comps_info %>% 
           filter(Class == "Alkane") %>% 
           pull(Chain.length)))
Ca_false_std_info

Ca_true_std_info <- std_info %>% 
  filter(!Chain.length %in% (Ca_true_comps_info %>% 
          filter(Class == "Alkane") %>% 
          pull(Chain.length)))
Ca_true_std_info

Ib_false_std_info <- std_info %>% 
  filter(!Chain.length %in% (Ib_false_comps_info %>% 
                               filter(Class == "Alkane") %>% 
                               pull(Chain.length)))
Ib_false_std_info

Ib_true_std_info <- std_info %>% 
  filter(!Chain.length %in% (Ib_true_comps_info %>% 
                               filter(Class == "Alkane") %>% 
                               pull(Chain.length)))
Ib_true_std_info

## Calculate RI ####
# Merge comps_info and std_info, so the mean_Rt of every alkane 
# is in the same data frame
Ca_false_comps_info <- merge(Ca_false_std_info
                       , Ca_false_comps_info
                       , all = T
                       , sort = F) %>%
  arrange(mean_RT) %>% 
  select(Peak
         , everything()) %>% 
  as_tibble()
Ca_false_comps_info

Ca_true_comps_info <- merge(Ca_true_std_info
                       , Ca_true_comps_info
                       , all = T
                       , sort = F) %>%
  arrange(mean_RT) %>% 
  select(Peak
         , everything()) %>% 
  as_tibble()
Ca_true_comps_info

Ib_false_comps_info <- merge(Ib_false_std_info
                            , Ib_false_comps_info
                            , all = T
                            , sort = F) %>%
  arrange(mean_RT) %>% 
  select(Peak
         , everything()) %>% 
  as_tibble()
Ib_false_comps_info

Ib_true_comps_info <- merge(Ib_true_std_info
                            , Ib_true_comps_info
                            , all = T
                            , sort = F) %>%
  arrange(mean_RT) %>% 
  select(Peak
         , everything()) %>% 
  as_tibble()
Ib_true_comps_info

# Use the k.ri() function to calculate the retention index
Ca_false_comps_info <- k.ri(Ca_false_comps_info)
Ca_true_comps_info <- k.ri(Ca_true_comps_info)
Ib_false_comps_info <- k.ri(Ib_false_comps_info)
Ib_true_comps_info <- k.ri(Ib_true_comps_info)

# Relative abundance  (%) ----
## Extract the abundance data into a data frame
Ca_false_daten <- Ca_false_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ca_true_daten <- Ca_true_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ib_false_daten <- Ib_false_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ib_true_daten <- Ib_true_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

## Replace NAs with 0s
Ca_false_daten[is.na(Ca_false_daten)] <- 0
Ca_true_daten[is.na(Ca_true_daten)] <- 0
Ib_false_daten[is.na(Ib_false_daten)] <- 0
Ib_true_daten[is.na(Ib_true_daten)] <- 0



### Calculate the relative abundance (%) of each compound per sample 
Ca_false_daten <- Ca_false_daten %>% t / rowSums(Ca_false_daten %>% t) * 100
Ca_false_daten <- Ca_false_daten %>% 
  t %>% 
  as.data.frame() %>% 
  as_tibble()

Ca_true_daten <- Ca_true_daten %>% t / rowSums(Ca_true_daten %>% t) * 100
Ca_true_daten <- Ca_true_daten %>% 
  t %>% 
  as.data.frame() %>% 
  as_tibble()

Ib_false_daten <- Ib_false_daten %>% t / rowSums(Ib_false_daten %>% t) * 100
Ib_false_daten <- Ib_false_daten %>% 
  t %>% 
  as.data.frame() %>% 
  as_tibble()

Ib_true_daten <- Ib_true_daten %>% t / rowSums(Ib_true_daten %>% t) * 100
Ib_true_daten <- Ib_true_daten %>% 
  t %>% 
  as.data.frame() %>% 
  as_tibble()

### Verify that the sum of all relative abundances per sample is exactly 100
Ca_false_daten %>% colSums()
Ca_true_daten %>% colSums()
Ib_false_daten %>% colSums()
Ib_true_daten %>% colSums()

# Export final data frames ----
## Sort individuals' data columns according to the order of the individuals
## in the grouping_info data frames
Ca_false_daten <- Ca_false_daten %>% select(all_of(grouping_info %>% 
                                         filter(Subspecies == "Ca" 
                                                & Flying == "FALSE") %>% 
                                         pull(Bee_number)))

Ca_true_daten <- Ca_true_daten %>% select(all_of(grouping_info %>% 
                                         filter(Subspecies == "Ca" 
                                                & Flying == "TRUE") %>% 
                                         pull(Bee_number)))

Ib_false_daten <- Ib_false_daten %>% select(all_of(grouping_info %>% 
                                                   filter(Subspecies == "Ib" 
                                                          & Flying == "FALSE") %>% 
                                                   pull(Bee_number)))

Ib_true_daten <- Ib_true_daten %>% select(all_of(grouping_info %>% 
                                                   filter(Subspecies == "Ib" 
                                                          & Flying == "TRUE") %>% 
                                                   pull(Bee_number)))

## Adjust the area data frame
Ca_false_area <- cbind(Ca_false_area %>% 
                   select(Peak:Compound)
                 , Ca_false_daten) %>% 
  as_tibble()

Ca_true_area <- cbind(Ca_true_area %>% 
                   select(Peak:Compound)
                 , Ca_true_daten) %>% 
  as_tibble()

Ib_false_area <- cbind(Ib_false_area %>% 
                        select(Peak:Compound)
                      , Ib_false_daten) %>% 
  as_tibble()

Ib_true_area <- cbind(Ib_true_area %>% 
                        select(Peak:Compound)
                      , Ib_true_daten) %>% 
  as_tibble()

## merge the area data frame and comps_info data frame to adjust peak numbering
Ca_false_table <- my_merge(Ca_false_comps_info
                  , Ca_false_area
                  , merge_by = c("Peak", "Compound")) %>% 
  as.data.frame()
Ca_false_table
str(Ca_false_table)

Ca_true_table <- my_merge(Ca_true_comps_info
                  , Ca_true_area
                  , merge_by = c("Peak", "Compound")) %>% 
  as.data.frame()
Ca_true_table
str(Ca_true_table)

Ib_false_table <- my_merge(Ib_false_comps_info
                          , Ib_false_area
                          , merge_by = c("Peak", "Compound")) %>% 
  as.data.frame()
Ib_false_table
str(Ib_false_table)

Ib_true_table <- my_merge(Ib_true_comps_info
                          , Ib_true_area
                          , merge_by = c("Peak", "Compound")) %>% 
  as.data.frame()
Ib_true_table
str(Ib_true_table)

## Adjust peak numbering
Ca_false_table$Peak <- paste0("P", 1:length(Ca_false_table$Peak))
head(Ca_false_table)
str(Ca_false_table)

Ca_true_table$Peak <- paste0("P", 1:length(Ca_true_table$Peak))
head(Ca_true_table)
str(Ca_true_table)

Ib_false_table$Peak <- paste0("P", 1:length(Ib_false_table$Peak))
head(Ib_false_table)
str(Ib_false_table)

Ib_true_table$Peak <- paste0("P", 1:length(Ib_true_table$Peak))
head(Ib_true_table)
str(Ib_true_table)

## Define peak numbers as row names
rownames(Ca_false_table) <- Ca_false_table$Peak
Ca_false_table <- Ca_false_table %>% select(-Peak)
head(Ca_false_table)
str(Ca_false_table)

rownames(Ca_true_table) <- Ca_true_table$Peak
Ca_true_table <- Ca_true_table %>% select(-Peak)
head(Ca_true_table)
str(Ca_true_table)

rownames(Ib_false_table) <- Ib_false_table$Peak
Ib_false_table <- Ib_false_table %>% select(-Peak)
head(Ib_false_table)
str(Ib_false_table)

rownames(Ib_true_table) <- Ib_true_table$Peak
Ib_true_table <- Ib_true_table %>% select(-Peak)
head(Ib_true_table)
str(Ib_true_table)

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
save(list = c("Ca_false_table"
              , "Ca_true_table"
              , "Ib_false_table"
              , "Ib_true_table"
              , "grouping_info")
     , file = here("data"
                   , "processed"
                   , "Flying"
                   , "group_tables.Rdata"))
print("The data frames for analysis were exported")

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "Flying", "SInf_Script04.txt"))
print("The sessionInfo report was exported. The script 04 finished running")

# ## Detach/unload packages
# lapply(NPacks, unloadNamespace)
# 
# ## Clear environment
# rm(list=ls())
