

# Install and/or load packages ----
NPacks <- c("tidyr", "dplyr", "here")

#install.packages(c(NPacks,"pacman"))

## Load all the packages and try to install them if they are not
## (Uses/require pacman package)
pacman::p_load(char = NPacks)

# Load the nestmate recognition data ----
nm.rec_daten <- read.csv(here("data", "raw", "NM-recognition"
                              , "nestmate_recognition-data.csv"))

# Re-structure NM recognition data set 1/2 ----
nm.rec_daten <- nm.rec_daten %>% 
  ## Add group variable indicating whether the bee was a nest mate
  ## or not of the guards at the hive were it was presented
  mutate(group = ifelse(nm.rec_daten$Source.hive == nm.rec_daten$Target.hive
                        , "nestmate"
                        , "non.nestmate")) %>% 
  ## Removes columns that indicates video file title and time of the day
  select(!Video & !Time) %>%
  ## Reorders the columns, leaving acceptance at the end
  select(Sample.ID
         , Date:Target.hive
         , group
         , Acceptance) %>%
  as.data.frame()

## This transforms the date into the year only
nm.rec_daten$Date <- nm.rec_daten$Date %>% as.Date() %>% format('%Y') 

nm.rec_daten$location[nm.rec_daten$location == "WUE"] <- "Würzburg"
nm.rec_daten$location[nm.rec_daten$location == "BRA"] <- "Bragança"

nm.rec_daten$Subspecies[nm.rec_daten$Subspecies == "Ib"] <- "A. m. iberiensis"
nm.rec_daten$Subspecies[nm.rec_daten$Subspecies == "Ca"] <- "A. m. carnica"

str(nm.rec_daten)

#  ----
nm.rec_daten %>% group_by(Subspecies) %>%
        summarise(mean = mean(Acceptance), n=n()) %>% 
  mutate(T = n*mean, F = n - n*mean) %>% 
  as.data.frame()

nm.rec_daten %>% filter(Subspecies == "A. m. iberiensis")  %>% group_by(group) %>% 
  summarise(mean = mean(Acceptance), n=n()) %>% 
  mutate(T = n*mean, F = n - n*mean) %>% 
  as.data.frame()

nm.rec_daten %>% filter(Subspecies == "A. m. carnica") %>% group_by(group) %>% 
  summarise(mean = mean(Acceptance), n=n()) %>% 
  mutate(T = n*mean, F = n - n*mean) %>% 
  as.data.frame()

# Acceptance data set ----
accept_daten <- nm.rec_daten %>% group_by(location, Subspecies
                                          , Target.hive, group) %>% 
  count(Acceptance) %>%
 # mutate(hive=paste(location, Subspecies, Target.hive, sep = "_")) %>%
  pivot_wider(id_cols = c(location, Subspecies, Target.hive)
              , names_from = c(group, Acceptance), values_from = n) %>% 
  as.data.frame()

accept_daten[is.na(accept_daten)] <- 0

accept_daten <- accept_daten %>% 
  mutate(n_nestmates = nestmate_TRUE + nestmate_FALSE
         ,n_non.nestmates = non.nestmate_TRUE + non.nestmate_FALSE) %>%
  as.data.frame()

accept_daten[, colnames(accept_daten %>% 
                          select(location:Target.hive))] <- 
  lapply(accept_daten[, colnames(accept_daten %>% 
                                   select(location:Target.hive))], as.factor)

accept_daten[, colnames(accept_daten %>% 
                          select(!where(is.factor)))] <- 
  lapply(accept_daten[, colnames(accept_daten %>% 
                                   select(!where(is.factor)))], as.integer)

colnames(accept_daten)[colnames(accept_daten) == "Target.hive"] <- "Hive"

accept_daten <- accept_daten %>% 
  mutate(nestmates = nestmate_TRUE/n_nestmates
         , non.nestmates = non.nestmate_TRUE/n_non.nestmates) %>% 
  #select(Hive, everything()) %>% select(Hive:Subspecies, n_nestmates:non.nestmates) %>% 
  select(-(nestmate_TRUE:nestmate_FALSE)) %>% 
  as.data.frame()

str(accept_daten)

# Re-structure NM recognition data set 2/2 ----
rownames(nm.rec_daten) <- nm.rec_daten$Sample.ID
nm.rec_daten <- nm.rec_daten %>% select(!Sample.ID) %>% as.data.frame()

nm.rec_daten[, colnames(nm.rec_daten %>% 
                         select(!"Date" & !"Acceptance"))] <- 
  lapply(nm.rec_daten[, colnames(nm.rec_daten %>% 
                                select(!"Date" & !"Acceptance"))], as.factor)

str(nm.rec_daten)

# Export data sets ----
write.csv(nm.rec_daten, here("data", "processed", "NM-Recognition"
                             , "nestmate_recognition-data_set.csv"))
saveRDS(nm.rec_daten, here("data", "processed", "NM-Recognition"
                           , "nm.rec_daten.Rds"))

write.csv(accept_daten, here("data", "processed", "NM-Recognition"
                             , "acceptance-data_set.csv"))
saveRDS(accept_daten, here("data", "processed", "NM-Recognition"
                           , "accept_daten.Rds"))

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output"
                             , "sInf_NM_Recognition-script_01.txt"))

## Detach/unload packages
lapply(NPacks, unloadNamespace)

