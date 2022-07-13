

# Install and/or load packages ----
NPacks <- c("here", "tidyr", "dplyr")

#install.packages(c(NPacks,"pacman"))

## Load all the packages and try to install them if they are not (Uses/require pacman package)
pacman::p_load(char=NPacks)

# Load data ----
s.data<-read.csv(here("data", "raw"#, "Survival"
                      , "Survival experiment_11_07_2022.csv"  #"Survival-data.csv"
                      )
                 ,row.names=1,check.names=FALSE)
 s.data <- s.data %>% filter(!is.na(Status)) #%>% select(!GCMS.file)
# s.data$Time <- as.numeric(s.data$Time)

s.data[, "Season"][s.data[, "Season"] == "S"] <- "Summer"
#s.data[, "Season"][s.data[, "Season"] == "W"] <- "Winter"  
s.data[, "Task"][s.data[, "Task"] == "Nu"] <- "Nurses"
#s.data[, "Task"][s.data[, "Task"] == "NPF"] <- "Non-pollen foragers"
s.data[, "Task"][s.data[, "Task"] == "PF"] <- "Pollen foragers"
s.data[, "Task"][s.data[, "Task"] == "IW"] <- "In-hive workers"
s.data[, "Task"][s.data[, "Task"] == "OW"] <- "Out-hive workers"
s.data[, "Subsp."][s.data[, "Subsp."] == "Ib"] <- "A. m. iberiensis"
s.data[, "Subsp."][s.data[, "Subsp."] == "Ca"] <- "A. m. carnica"
#s.data[, "Temperature"][s.data[, "Temperature"] == ""] <- ""
#s.data[, "Temperature"][s.data[, "Temperature"] == ""] <- ""

s.data <- s.data %>%
  mutate(Season = ifelse(s.data$Season == "S"
                         , "Summer"
                         , "Winter")
         , Temperature = ifelse(s.data$Temp. == "45"
                                , "45°C"
                                , "20°C")
         , Silica = ifelse(s.data$Silica == T
                           , "with silica"
                           , "without silica"))
  
s.data[, colnames(s.data %>% 
                    select(Season:Silica, Temperature))] <- 
  lapply(s.data[, colnames(s.data %>% 
                             select(Season:Silica, Temperature))], as.factor)
str(s.data)

# Export data sets ----
#write.csv(s.data, "data/processed/Survival/Survival-data.csv")
saveRDS(s.data, here("data", "raw"#, "Survival"
                     , "Survival-data.Rds"))

# End ----
## Report session information
capture.output(sessionInfo(), file = "output/sInf_Survival-script_01.txt")

## Detach/unload packages
lapply(NPacks, unloadNamespace)
