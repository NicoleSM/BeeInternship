
# Install and/or load packages ----
NPacks <- c("GGally", "here", "tidyverse")

#install.packages(c(NPacks,"pacman"))

## Load all the packages and try to install them if they are not (Uses/require pacman package)
pacman::p_load(char=NPacks)

# Load data ----
s.data <- readRDS(here("data", "processed", "Survival", "survival-data.Rds"))

str(s.data)

s.data %>% 
  select(-all_of(c("Date"
                   , "Collection.Hive"
                   , "Native.Hive"
                   , "Fostered"
                   , "Status"))) %>% 
  select(Season:Subspecies) %>% 
  ggpairs(aes(color = '#')) +
  scale_fill_viridis_d() +
  #coord_flip() +
  theme_minimal()

s.data %>% filter(Season == "Summer") %>% 
  select(-all_of(c("Date"
                   , "Collection.Hive"
                   , "Native.Hive"
                   , "Fostered"
                   , "Status"))) %>% 
  ggpairs(aes(color = Task
              , alpha = 0.95)) +
  scale_fill_viridis_d()

s.data %>% filter(Season == "Winter" & !Collection.Hive==3) %>% 
  select(-all_of(c("Date"
                   , "Collection.Hive"
                   , "Native.Hive"
                   , "Fostered"
                   , "Status"))) %>% 
  ggpairs(aes(color = Task
              , alpha = 0.95)) +
  scale_fill_viridis_d()

