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
## Load aligned GC-MS tables generated by script 03
load(here("data", "processed", "Flying", "aligned_gcms-data.Rdata"))

## Load the compound information, result from the identification process
path_comps_id <- list.files(path = here("data", "raw", "Flying", "tmp")
                            , pattern = "ID_table.csv"
                            , full.names = T)
comps_id <- lapply(path_comps_id, read_csv)
summary(comps_id)

names(comps_id) <- str_split(path_comps_id
                                 , "/"
                                 , simplify = T) %>% 
  str_subset(".csv") %>% 
  str_remove("ID_table.csv") 
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
### Ca-false ####
### Extract compounds information from comps_id
Ca-false_comps_info <- cbind(comps_id$Ca-false 
                  , data.frame(row.names = rownames(comps_id$Ca-false)
                               , t(as.data.frame(strsplit(comps_id$Ca-false$Compound
                                                          , "_"))))) %>% 
  as_tibble()

### Define columns names
colnames(Ca-false_comps_info) <- c(colnames(comps_id$Ca-false)
                          , "Chain.length", "Class", "Mod.position")
Ca-false_comps_info

### Correct entries format
#### Compound names
# Store the compound names in a vector, where they will be altered into their final format
comps_names <- Ca-false_comps_info$Compound

# Set the first part of the identified CHC names (full name in the case of the alkanes)
comps_names[!is.na(comps_names)] <- paste0(Ca-false_comps_info %>% 
                                             filter(!is.na(Compound)) %>% 
                                             pull(Mod.position)
                                           , Ca-false_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Class)
                                           , Ca-false_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Chain.length)) %>% 
  str_replace_all(paste(c("ane", "ene", "diene"), collapse = '|'), "NA") %>% 
  str_remove_all("NA")

# Finish formatting the names of the identified CHC
comps_names[!is.na(comps_names)] <- paste(comps_names[!is.na(comps_names)]
                                          , Ca-false_comps_info %>% 
                                            filter(!is.na(Compound)) %>% 
                                            pull(Class) %>% 
                                            str_replace(paste(c("ane", "Me"
                                                                , "Dime", "Trime")
                                                              , collapse = '|')
                                                        , "NA")
                                          , sep = ":") %>% 
  str_remove(":NA")

# Change the compound names in the comps_info data frame to the correct format names
Ca-false_comps_info$Compound <- comps_names
Ca-false_comps_info
Ca-false_comps_info %>% filter(!is.na(Compound))

#### Chain length
Ca-false_comps_info$Chain.length <- Ca-false_comps_info$Chain.length %>%
  str_remove("C") %>%
  as.integer()
Ca-false_comps_info
Ca-false_comps_info %>% filter(!is.na(Compound))

#### Class
Ca-false_comps_info['Class'][Ca-false_comps_info['Class'] == "ane"] <- "Alkane"
Ca-false_comps_info['Class'][Ca-false_comps_info['Class'] == "ene"] <- "Alkene"
Ca-false_comps_info['Class'][Ca-false_comps_info['Class'] == "diene"] <- "Alkadiene"
Ca-false_comps_info['Class'][Ca-false_comps_info['Class'] == "Me"] <- "Methyl"
Ca-false_comps_info['Class'][Ca-false_comps_info['Class'] == "Dime"] <- "Dimethyl"
Ca-false_comps_info['Class'][Ca-false_comps_info['Class'] == "Trime"] <- "Trimethyl"
Ca-false_comps_info['Class'][Ca-false_comps_info['Class'] == "Tetrame"] <- "Tetramethyl"
Ca-false_comps_info
Ca-false_comps_info %>% filter(!is.na(Compound))

#### Add mean_RT 
Ca-false_comps_info <- my_merge(Ca-false_comps_info
                       , Ca-false_RT %>% 
                         select(Peak, mean_RT)
                       , merge_by = "Peak") %>% 
  arrange(mean_RT) %>% 
  as_tibble()
Ca-false_comps_info
Ca-false_comps_info %>% filter(!is.na(Compound))

### Ca-true ####
### Extract compounds information from comps_id
Ca-true_comps_info <- cbind(comps_id$Ca-true
                       , data.frame(row.names = rownames(comps_id$Ca-true)
                                    , t(as.data.frame(strsplit(comps_id$Ca-true$Compound
                                                               , "_"))))) %>% 
  as_tibble()

### Define columns names
colnames(Ca-true_comps_info) <- c(colnames(comps_id$Ca-true)
                             , "Chain.length", "Class", "Mod.position")
Ca-true_comps_info

### Correct entries format
#### Compound names
# Store the compound names in a vector, where they will be altered into their final format
comps_names <- Ca-true_comps_info$Compound

# Set the first part of the identified CHC names (full name in the case of the alkanes)
comps_names[!is.na(comps_names)] <- paste0(Ca-true_comps_info %>% 
                                             filter(!is.na(Compound)) %>% 
                                             pull(Mod.position)
                                           , Ca-true_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Class)
                                           , Ca-true_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Chain.length)) %>% 
  str_replace_all(paste(c("ane", "ene", "diene"), collapse = '|'), "NA") %>% 
  str_remove_all("NA")

# Finish formatting the names of the identified CHC
comps_names[!is.na(comps_names)] <- paste(comps_names[!is.na(comps_names)]
                                          , Ca-true_comps_info %>% 
                                            filter(!is.na(Compound)) %>% 
                                            pull(Class) %>% 
                                            str_replace(paste(c("ane", "Me")
                                                              , collapse = '|')
                                                        , "NA")
                                          , sep = ":") %>% 
  str_remove(":NA")

# Change the compound names in the comps_info data frame to the correct format names
Ca-true_comps_info$Compound <- comps_names
Ca-true_comps_info
Ca-true_comps_info %>% filter(!is.na(Compound))

#### Chain length
Ca-true_comps_info$Chain.length <- Ca-true_comps_info$Chain.length %>%
  str_remove("C") %>%
  as.integer()
Ca-true_comps_info
Ca-true_comps_info %>% filter(!is.na(Compound))

#### Class
Ca-true_comps_info['Class'][Ca-true_comps_info['Class'] == "ane"] <- "Alkane"
Ca-true_comps_info['Class'][Ca-true_comps_info['Class'] == "ene"] <- "Alkene"
Ca-true_comps_info['Class'][Ca-true_comps_info['Class'] == "diene"] <- "Alkadiene"
Ca-true_comps_info['Class'][Ca-true_comps_info['Class'] == "Me"] <- "Methyl"
Ca-true_comps_info['Class'][Ca-true_comps_info['Class'] == "Dime"] <- "Dimethyl"
Ca-true_comps_info['Class'][Ca-true_comps_info['Class'] == "Trime"] <- "Trimethyl"
Ca-true_comps_info['Class'][Ca-true_comps_info['Class'] == "Tetrame"] <- "Tetramethyl"
Ca-true_comps_info
Ca-true_comps_info %>% filter(!is.na(Compound))

#### Add mean_RT 
Ca-true_comps_info <- my_merge(Ca-true_comps_info
                          , Ca-true_RT %>% 
                            select(Peak, mean_RT)
                          , merge_by = "Peak") %>% 
  arrange(mean_RT) %>% 
  as_tibble()
Ca-true_comps_info
Ca-true_comps_info %>% filter(!is.na(Compound))

### Ib-false ####
### Extract compounds information from comps_id
Ib-false_comps_info <- cbind(comps_id$Ib-false
                            , data.frame(row.names = rownames(comps_id$Ib-false)
                                         , t(as.data.frame(strsplit(comps_id$Ib-false$Compound
                                                                    , "_"))))) %>% 
  as_tibble()

### Define columns names
colnames(Ib-false_comps_info) <- c(colnames(comps_id$Ib-false)
                                  , "Chain.length", "Class", "Mod.position")
Ib-false_comps_info

### Correct entries format
#### Compound names
# Store the compound names in a vector, where they will be altered into their final format
comps_names <- Ib-false_comps_info$Compound

# Set the first part of the identified CHC names (full name in the case of the alkanes)
comps_names[!is.na(comps_names)] <- paste0(Ib-false_comps_info %>% 
                                             filter(!is.na(Compound)) %>% 
                                             pull(Mod.position)
                                           , Ib-false_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Class)
                                           , Ib-false_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Chain.length)) %>% 
  str_replace_all(paste(c("ane", "ene", "diene"), collapse = '|'), "NA") %>% 
  str_remove_all("NA")

# Finish formatting the names of the identified CHC
comps_names[!is.na(comps_names)] <- paste(comps_names[!is.na(comps_names)]
                                          , Ib-false_comps_info %>% 
                                            filter(!is.na(Compound)) %>% 
                                            pull(Class) %>% 
                                            str_replace(paste(c("ane", "Me")
                                                              , collapse = '|')
                                                        , "NA")
                                          , sep = ":") %>% 
  str_remove(":NA")

# Change the compound names in the comps_info data frame to the correct format names
Ib-false_comps_info$Compound <- comps_names
Ib-false_comps_info
Ib-false_comps_info %>% filter(!is.na(Compound))

#### Chain length
Ib-false_comps_info$Chain.length <- Ib-false_comps_info$Chain.length %>%
  str_remove("C") %>%
  as.integer()
Ib-false_comps_info
Ib-false_comps_info %>% filter(!is.na(Compound))

#### Class
Ib-false_comps_info['Class'][Ib-false_comps_info['Class'] == "ane"] <- "Alkane"
Ib-false_comps_info['Class'][Ib-false_comps_info['Class'] == "ene"] <- "Alkene"
Ib-false_comps_info['Class'][Ib-false_comps_info['Class'] == "diene"] <- "Alkadiene"
Ib-false_comps_info['Class'][Ib-false_comps_info['Class'] == "Me"] <- "Methyl"
Ib-false_comps_info['Class'][Ib-false_comps_info['Class'] == "Dime"] <- "Dimethyl"
Ib-false_comps_info['Class'][Ib-false_comps_info['Class'] == "Trime"] <- "Trimethyl"
Ib-false_comps_info['Class'][Ib-false_comps_info['Class'] == "Tetrame"] <- "Tetramethyl"
Ib-false_comps_info
Ib-false_comps_info %>% filter(!is.na(Compound))

#### Add mean_RT 
Ib-false_comps_info <- my_merge(Ib-false_comps_info
                               , Ib-false_RT %>% 
                                 select(Peak, mean_RT)
                               , merge_by = "Peak") %>% 
  arrange(mean_RT) %>% 
  as_tibble()
Ib-false_comps_info
Ib-false_comps_info %>% filter(!is.na(Compound))

### Ib-true ####
### Extract compounds information from comps_id
Ib-true_comps_info <- cbind(comps_id$Ib-true
                            , data.frame(row.names = rownames(comps_id$Ib-true)
                                         , t(as.data.frame(strsplit(comps_id$Ib-true$Compound
                                                                    , "_"))))) %>% 
  as_tibble()

### Define columns names
colnames(Ib-true_comps_info) <- c(colnames(comps_id$Ib-true)
                                  , "Chain.length", "Class", "Mod.position")
Ib-true_comps_info

### Correct entries format
#### Compound names
# Store the compound names in a vector, where they will be altered into their final format
comps_names <- Ib-true_comps_info$Compound

# Set the first part of the identified CHC names (full name in the case of the alkanes)
comps_names[!is.na(comps_names)] <- paste0(Ib-true_comps_info %>% 
                                             filter(!is.na(Compound)) %>% 
                                             pull(Mod.position)
                                           , Ib-true_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Class)
                                           , Ib-true_comps_info %>%
                                             filter(!is.na(Compound)) %>% 
                                             pull(Chain.length)) %>% 
  str_replace_all(paste(c("ane", "ene", "diene"), collapse = '|'), "NA") %>% 
  str_remove_all("NA")

# Finish formatting the names of the identified CHC
comps_names[!is.na(comps_names)] <- paste(comps_names[!is.na(comps_names)]
                                          , Ib-true_comps_info %>% 
                                            filter(!is.na(Compound)) %>% 
                                            pull(Class) %>% 
                                            str_replace(paste(c("ane", "Me")
                                                              , collapse = '|')
                                                        , "NA")
                                          , sep = ":") %>% 
  str_remove(":NA")

# Change the compound names in the comps_info data frame to the correct format names
Ib-true_comps_info$Compound <- comps_names
Ib-true_comps_info
Ib-true_comps_info %>% filter(!is.na(Compound))

#### Chain length
Ib-true_comps_info$Chain.length <- Ib-true_comps_info$Chain.length %>%
  str_remove("C") %>%
  as.integer()
Ib-true_comps_info
Ib-true_comps_info %>% filter(!is.na(Compound))

#### Class
Ib-true_comps_info['Class'][Ib-true_comps_info['Class'] == "ane"] <- "Alkane"
Ib-true_comps_info['Class'][Ib-true_comps_info['Class'] == "ene"] <- "Alkene"
Ib-true_comps_info['Class'][Ib-true_comps_info['Class'] == "diene"] <- "Alkadiene"
Ib-true_comps_info['Class'][Ib-true_comps_info['Class'] == "Me"] <- "Methyl"
Ib-true_comps_info['Class'][Ib-true_comps_info['Class'] == "Dime"] <- "Dimethyl"
Ib-true_comps_info['Class'][Ib-true_comps_info['Class'] == "Trime"] <- "Trimethyl"
Ib-true_comps_info['Class'][Ib-true_comps_info['Class'] == "Tetrame"] <- "Tetramethyl"
Ib-true_comps_info
Ib-true_comps_info %>% filter(!is.na(Compound))

#### Add mean_RT 
Ib-true_comps_info <- my_merge(Ib-true_comps_info
                               , Ib-true_RT %>% 
                                 select(Peak, mean_RT)
                               , merge_by = "Peak") %>% 
  arrange(mean_RT) %>% 
  as_tibble()
Ib-true_comps_info
Ib-true_comps_info %>% filter(!is.na(Compound))

## Samples ####
### Ca-false ####
### Add CHC identification to the aligned sample data frames
Ca-false_RT <- my_merge(Ca-false_comps_info %>% 
                    select(Peak, Compound)
                  , Ca-false_RT %>% 
                    select(!mean_RT)
                  , merge_by = "Peak") %>%
  as_tibble()
Ca-false_RT

Ca-false_area <- my_merge(Ca-false_comps_info %>% 
                      select(Peak, Compound)
                    , Ca-false_area %>% 
                      select(!mean_RT)
                    , merge_by = "Peak") %>%
  as_tibble()
Ca-false_area

### Ca-true ####
### Add CHC identification to the aligned sample data frames
Ca-true_RT <- my_merge(Ca-true_comps_info %>% 
                    select(Peak, Compound)
                  , Ca-true_RT %>% 
                    select(!mean_RT)
                  , merge_by = "Peak") %>%
  as_tibble()
Ca-true_RT

Ca-true_area <- my_merge(Ca-true_comps_info %>% 
                      select(Peak, Compound)
                    , Ca-true_area %>% 
                      select(!mean_RT)
                    , merge_by = "Peak") %>%
  as_tibble()
Ca-true_area

### Ib-false ####
### Add CHC identification to the aligned sample data frames
Ib-false_RT <- my_merge(Ib-false_comps_info %>% 
                         select(Peak, Compound)
                       , Ib-false_RT %>% 
                         select(!mean_RT)
                       , merge_by = "Peak") %>%
  as_tibble()
Ib-false_RT

Ib-false_area <- my_merge(Ib-false_comps_info %>% 
                           select(Peak, Compound)
                         , Ib-false_area %>% 
                           select(!mean_RT)
                         , merge_by = "Peak") %>%
  as_tibble()
Ib-false_area

### Ib-true ####
### Add CHC identification to the aligned sample data frames
Ib-true_RT <- my_merge(Ib-true_comps_info %>% 
                         select(Peak, Compound)
                       , Ib-true_RT %>% 
                         select(!mean_RT)
                       , merge_by = "Peak") %>%
  as_tibble()
Ib-true_RT

Ib-true_area <- my_merge(Ib-true_comps_info %>% 
                           select(Peak, Compound)
                         , Ib-true_area %>% 
                           select(!mean_RT)
                         , merge_by = "Peak") %>%
  as_tibble()
Ib-true_area

### grouping_info ####
### Define grouping variables as factors
grouping_info[, colnames(grouping_info)] <-
  lapply(grouping_info[, colnames(grouping_info)], as.factor)
grouping_info

# Clean the data ----
## Delete trace compounds ####
### Extract the abundance data into a data frame
Ca-false_daten <- Ca-false_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ca-true_daten <- Ca-true_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ib-false_daten <- Ib-false_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ib-true_daten <- Ib-true_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

### Use trace_comps function to delete compounds, from each sample, 
### with an abundance below a defined threshold
Ca-false_daten <- trace_comps(Ca-false_daten, threshold = 0.01)
Ca-true_daten <- trace_comps(Ca-true_daten, threshold = 0.01)
Ib-false_daten <- trace_comps(Ib-false_daten, threshold = 0.01)
Ib-true_daten <- trace_comps(Ib-true_daten, threshold = 0.01)

### Adjust the area data frame
Ca-false_area <- cbind(Ca-false_area %>% 
                   select(Peak:Compound)
                 , Ca-false_daten) %>% 
  as_tibble()
Ca-false_area

Ca-true_area <- cbind(Ca-true_area %>% 
                   select(Peak:Compound)
                 , Ca-true_daten) %>% 
  as_tibble()
Ca-true_area

Ib-false_area <- cbind(Ib-false_area %>% 
                        select(Peak:Compound)
                      , Ib-false_daten) %>% 
  as_tibble()
Ib-false_area

Ib-true_area <- cbind(Ib-true_area %>% 
                        select(Peak:Compound)
                      , Ib-true_daten) %>% 
  as_tibble()
Ib-true_area

### Adjust the RT data frame accordingly
Ca-false_RT_daten <- Ca-false_RT %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ca-true_RT_daten <- Ca-true_RT %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ib-false_RT_daten <- Ib-false_RT %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ib-true_RT_daten <- Ib-true_RT %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ca-false_RT_daten[is.na(Ca-false_daten)] <- NA
Ca-false_RT <- cbind(Ca-false_RT %>% 
                 select(Peak:Compound)
               , Ca-false_RT_daten) %>% 
  as_tibble()
Ca-false_RT

Ca-true_RT_daten[is.na(Ca-true_daten)] <- NA
Ca-true_RT <- cbind(Ca-true_RT %>% 
                 select(Peak:Compound)
               , Ca-true_RT_daten) %>% 
  as_tibble()
Ca-true_RT

Ib-false_RT_daten[is.na(Ib-false_daten)] <- NA
Ib-false_RT <- cbind(Ib-false_RT %>% 
                      select(Peak:Compound)
                    , Ib-false_RT_daten) %>% 
  as_tibble()
Ib-false_RT

Ib-true_RT_daten[is.na(Ib-true_daten)] <- NA
Ib-true_RT <- cbind(Ib-true_RT %>% 
                      select(Peak:Compound)
                    , Ib-true_RT_daten) %>% 
  as_tibble()
Ib-true_RT

## Delete non-CHC compounds ####
### Samples 
Ca-false_area <- Ca-false_area %>% filter(!is.na(Compound))
Ca-false_area
Ca-false_RT <- Ca-false_RT %>% filter(!is.na(Compound))
Ca-false_RT

Ca-true_area <- Ca-true_area %>% filter(!is.na(Compound))
Ca-true_area
Ca-true_RT <- Ca-true_RT %>% filter(!is.na(Compound))
Ca-true_RT

Ib-false_RT_daten[is.na(Ib-false_daten)] <- NA
Ib-false_RT <- cbind(Ib-false_RT %>% 
                      select(Peak:Compound)
                    , Ib-false_RT_daten) %>% 
  as_tibble()
Ib-false_RT

Ib-true_RT_daten[is.na(Ib-true_daten)] <- NA
Ib-true_RT <- cbind(Ib-true_RT %>% 
                      select(Peak:Compound)
                    , Ib-true_RT_daten) %>% 
  as_tibble()
Ib-true_RT

### Compounds information
Ca-false_comps_info <- Ca-false_comps_info %>% filter(!is.na(Compound))
Ca-false_comps_info

Ca-true_comps_info <- Ca-true_comps_info %>% filter(!is.na(Compound))
Ca-true_comps_info

Ib-true_comps_info <- Ib-true_comps_info %>% filter(!is.na(Compound))
Ib-true_comps_info

Ib-false_comps_info <- Ib-false_comps_info %>% filter(!is.na(Compound))
Ib-false_comps_info

### Standards
std_info <- std_info %>% filter(!is.na(Compound))
std_info

## Delete rare compounds - OPTIONAL ####
### Ca-false ####
Ca-false_area <- below_50(Ca-false_area
                    , grouping.info = grouping_info %>% 
                      filter(Task == "Non-flying Carnica pollen foragers"))
Ca-false_area

Ca-false_RT <- Ca-false_RT %>% 
  filter(Peak %in% Ca-false_area$Peak)
Ca-false_RT

Ca-false_comps_info <- Ca-false_comps_info %>% 
  filter(Peak %in% Ca-false_area$Peak)
Ca-false_comps_info

### Ca-true ####
Ca-true_area <- belCa-true_50(Ca-true_area
                    , grouping.info = grouping_info %>% 
                      filter(Task == "Flying Carnica pollen foragers"))
Ca-true_area

Ca-true_RT <- Ca-true_RT %>% 
  filter(Peak %in% Ca-true_area$Peak)
Ca-true_RT

Ca-true_comps_info <- Ca-true_comps_info %>% 
  filter(Peak %in% Ca-true_area$Peak)
Ca-true_comps_info

### Ib-false ####
Ib-false_area <- belIb-false_50(Ib-false_area
                              , grouping.info = grouping_info %>% 
                                filter(Task == "Non-flying Iberiensis pollen foragers"))
Ib-false_area

Ib-false_RT <- Ib-false_RT %>% 
  filter(Peak %in% Ib-false_area$Peak)
Ib-false_RT

Ib-false_comps_info <- Ib-false_comps_info %>% 
  filter(Peak %in% Ib-false_area$Peak)
Ib-false_comps_info

### Ib-true ####
Ib-true_area <- belIb-true_50(Ib-true_area
                              , grouping.info = grouping_info %>% 
                                filter(Task == "Flying Iberiensis pollen foragers"))
Ib-true_area

Ib-true_RT <- Ib-true_RT %>% 
  filter(Peak %in% Ib-true_area$Peak)
Ib-true_RT

Ib-true_comps_info <- Ib-true_comps_info %>% 
  filter(Peak %in% Ib-true_area$Peak)
Ib-true_comps_info

# Kováts Retention index ----
## Recalculate mean_RT ####
# It may have changed due to the fusion and deletion of peaks
# across the data set  
### Ca-false ####
identical(Ca-false_RT %>%
            select(!Peak:Compound) %>%
            rowMeans(na.rm = T)
          , Ca-false_comps_info$mean_RT)

Ca-false_comps_info$mean_RT <- Ca-false_RT %>%
  select(!Peak:Compound) %>% 
  rowMeans(na.rm = T)

identical(Ca-false_RT %>%
            select(!Peak:Compound) %>% 
            rowMeans(na.rm = T)
          , Ca-false_comps_info$mean_RT)
Ca-false_comps_info$mean_RT
Ca-false_comps_info

### Ca-true ####
identical(Ca-true_RT %>%
            select(!Peak:Compound) %>%
            rowMeans(na.rm = T)
          , Ca-true_comps_info$mean_RT)

Ca-true_comps_info$mean_RT <- Ca-true_RT %>%
  select(!Peak:Compound) %>% 
  rowMeans(na.rm = T)

identical(Ca-true_RT %>%
            select(!Peak:Compound) %>% 
            rowMeans(na.rm = T)
          , Ca-true_comps_info$mean_RT)
Ca-true_comps_info$mean_RT
Ca-true_comps_info

### Ib-false ####
identical(Ib-false_RT %>%
            select(!Peak:Compound) %>%
            rowMeans(na.rm = T)
          , Ib-false_comps_info$mean_RT)

Ib-false_comps_info$mean_RT <- Ib-false_RT %>%
  select(!Peak:Compound) %>% 
  rowMeans(na.rm = T)

identical(Ib-false_RT %>%
            select(!Peak:Compound) %>% 
            rowMeans(na.rm = T)
          , Ib-false_comps_info$mean_RT)
Ib-false_comps_info$mean_RT
Ib-false_comps_info


### Ib-true ####
identical(Ib-true_RT %>%
            select(!Peak:Compound) %>%
            rowMeans(na.rm = T)
          , Ib-true_comps_info$mean_RT)

Ib-true_comps_info$mean_RT <- Ib-true_RT %>%
  select(!Peak:Compound) %>% 
  rowMeans(na.rm = T)

identical(Ib-true_RT %>%
            select(!Peak:Compound) %>% 
            rowMeans(na.rm = T)
          , Ib-true_comps_info$mean_RT)
Ib-true_comps_info$mean_RT
Ib-true_comps_info

## Delete unnecessary standards ####
Ca-false_std_info <- std_info %>% 
  filter(!Chain.length %in% (Ca-false_comps_info %>% 
           filter(Class == "Alkane") %>% 
           pull(Chain.length)))
Ca-false_std_info

Ca-true_std_info <- std_info %>% 
  filter(!Chain.length %in% (Ca-true_comps_info %>% 
                               filter(Class == "Alkane") %>% 
                               pull(Chain.length)))
Ca-true_std_info

Ib-false_std_info <- std_info %>% 
  filter(!Chain.length %in% (Ib-false_comps_info %>% 
                               filter(Class == "Alkane") %>% 
                               pull(Chain.length)))
Ib-false_std_info

Ib-true_std_info <- std_info %>% 
  filter(!Chain.length %in% (Ib-true_comps_info %>% 
                               filter(Class == "Alkane") %>% 
                               pull(Chain.length)))
Ib-true_std_info

## Calculate RI ####
# Merge comps_info and std_info, so the mean_Rt of every alkane 
# is in the same data frame
Ca-false_comps_info <- merge(Ca-false_std_info
                       , Ca-false_comps_info
                       , all = T
                       , sort = F) %>%
  arrange(mean_RT) %>% 
  select(Peak
         , everything()) %>% 
  as_tibble()
Ca-false_comps_info

Ca-true_comps_info <- merge(Ca-true_std_info
                       , Ca-true_comps_info
                       , all = T
                       , sort = F) %>%
  arrange(mean_RT) %>% 
  select(Peak
         , everything()) %>% 
  as_tibble()
Ca-true_comps_info

Ib-false_comps_info <- merge(Ib-false_std_info
                            , Ib-false_comps_info
                            , all = T
                            , sort = F) %>%
  arrange(mean_RT) %>% 
  select(Peak
         , everything()) %>% 
  as_tibble()
Ib-false_comps_info

Ib-true_comps_info <- merge(Ib-true_std_info
                            , Ib-true_comps_info
                            , all = T
                            , sort = F) %>%
  arrange(mean_RT) %>% 
  select(Peak
         , everything()) %>% 
  as_tibble()
Ib-true_comps_info

# Use the k.ri() function to calculate the retention index
Ca-false_comps_info <- k.ri(Ca-false_comps_info)
Ca-true_comps_info <- k.ri(Ca-true_comps_info)
Ib-false_comps_info <- k.ri(Ib-false_comps_info)
Ib-true_comps_info <- k.ri(Ib-true_comps_info)

# Relative abundance  (%) ----
## Extract the abundance data into a data frame
Ca-false_daten <- Ca-false_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ca-true_daten <- Ca-true_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ib-false_daten <- Ib-false_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

Ib-true_daten <- Ib-true_area %>% 
  select(!Peak:Compound) %>%  
  as.data.frame()

## Replace NAs with 0s
Ca-false_daten[is.na(Ca-false_daten)] <- 0
Ca-true_daten[is.na(Ca-true_daten)] <- 0
Ib-false_daten[is.na(Ib-false_daten)] <- 0
Ib-true_daten[is.na(Ib-true_daten)] <- 0

### Calculate the relative abundance (%) of each compound per sample 
Ca-false_daten <- Ca-false_daten %>% t / rowSums(Ca-false_daten %>% t) * 100
Ca-false_daten <- Ca-false_daten %>% 
  t %>% 
  as.data.frame() %>% 
  as_tibble()

Ca-true_daten <- Ca-true_daten %>% t / rowSums(Ca-true_daten %>% t) * 100
Ca-true_daten <- Ca-true_daten %>% 
  t %>% 
  as.data.frame() %>% 
  as_tibble()

Ib-false_daten <- Ib-false_daten %>% t / rowSums(Ib-false_daten %>% t) * 100
Ib-false_daten <- Ib-false_daten %>% 
  t %>% 
  as.data.frame() %>% 
  as_tibble()

Ib-true_daten <- Ib-true_daten %>% t / rowSums(Ib-true_daten %>% t) * 100
Ib-true_daten <- Ib-true_daten %>% 
  t %>% 
  as.data.frame() %>% 
  as_tibble()

### Verify that the sum of all relative abundances per sample is exactly 100
Ca-false_daten %>% colSums()
Ca-true_daten %>% colSums()
Ib-false_daten %>% colSums()
Ib-true_daten %>% colSums()

# Export final data frames ----
## Sort individuals' data columns according to the order of the individuals
## in the grouping_info data frames
Ca-false_daten <- Ca-false_daten %>% select(all_of(grouping_info %>% 
                                         filter(Task == "Non-flying Carnica pollen foragers") %>% 
                                         pull(Individual)))

Ca-true_daten <- Ca-true_daten %>% select(all_of(grouping_info %>% 
                                         filter(Task == "Flying Carnica pollen foragers") %>% 
                                         pull(Individual)))

Ib-false_daten <- Ib-false_daten %>% select(all_of(grouping_info %>% 
                                                   filter(Task == "Non-flying Iberiensis pollen foragers") %>% 
                                                   pull(Individual)))

Ib-true_daten <- Ib-true_daten %>% select(all_of(grouping_info %>% 
                                                   filter(Task == "Flying Iberiensis pollen foragers") %>% 
                                                   pull(Individual)))

## Adjust the area data frame
Ca-false_area <- cbind(Ca-false_area %>% 
                   select(Peak:Compound)
                 , Ca-false_daten) %>% 
  as_tibble()

Ca-true_area <- cbind(Ca-true_area %>% 
                   select(Peak:Compound)
                 , Ca-true_daten) %>% 
  as_tibble()

Ib-false_area <- cbind(Ib-false_area %>% 
                        select(Peak:Compound)
                      , Ib-false_daten) %>% 
  as_tibble()

Ib-true_area <- cbind(Ib-true_area %>% 
                        select(Peak:Compound)
                      , Ib-true_daten) %>% 
  as_tibble()

## merge the area data frame and comps_info data frame to adjust peak numbering
Ca-false_table <- my_merge(Ca-false_comps_info
                  , Ca-false_area
                  , merge_by = c("Peak", "Compound")) %>% 
  as.data.frame()
Ca-false_table
str(Ca-false_table)

Ca-true_table <- my_merge(Ca-true_comps_info
                  , Ca-true_area
                  , merge_by = c("Peak", "Compound")) %>% 
  as.data.frame()
Ca-true_table
str(Ca-true_table)

Ib-false_table <- my_merge(Ib-false_comps_info
                          , Ib-false_area
                          , merge_by = c("Peak", "Compound")) %>% 
  as.data.frame()
Ib-false_table
str(Ib-false_table)

Ib-true_table <- my_merge(Ib-true_comps_info
                          , Ib-true_area
                          , merge_by = c("Peak", "Compound")) %>% 
  as.data.frame()
Ib-true_table
str(Ib-true_table)

## Adjust peak numbering
Ca-false_table$Peak <- paste0("P", 1:length(Ca-false_table$Peak))
head(Ca-false_table)
str(Ca-false_table)

Ca-true_table$Peak <- paste0("P", 1:length(Ca-true_table$Peak))
head(Ca-true_table)
str(Ca-true_table)

Ib-false_table$Peak <- paste0("P", 1:length(Ib-false_table$Peak))
head(Ib-false_table)
str(Ib-false_table)

Ib-true_table$Peak <- paste0("P", 1:length(Ib-true_table$Peak))
head(Ib-true_table)
str(Ib-true_table)

## Define peak numbers as row names
rownames(Ca-false_table) <- Ca-false_table$Peak
Ca-false_table <- Ca-false_table %>% select(-Peak)
head(Ca-false_table)
str(Ca-false_table)

rownames(Ca-true_table) <- Ca-true_table$Peak
Ca-true_table <- Ca-true_table %>% select(-Peak)
head(Ca-true_table)
str(Ca-true_table)

rownames(Ib-false_table) <- Ib-false_table$Peak
Ib-false_table <- Ib-false_table %>% select(-Peak)
head(Ib-false_table)
str(Ib-false_table)

rownames(Ib-true_table) <- Ib-true_table$Peak
Ib-true_table <- Ib-true_table %>% select(-Peak)
head(Ib-true_table)
str(Ib-true_table)

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
save(list = c("Ca-false_table"
              , "Ca-true_table"
              , "Ib-false_table"
              , "Ib-true_table"
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
print("The sessionInfo report was exported. The script 03 finished running")

# ## Detach/unload packages
# lapply(NPacks, unloadNamespace)
# 
# ## Clear environment
# rm(list=ls())
