
# Install and/or load packages ----
NPacks <- c("here", "tidyverse")

#install.packages(c(NPacks,"pacman"))

## Load all the packages and try to install them if they are not (Uses/require pacman package)
pacman::p_load(char = NPacks)

# Load data ----
s.data <- readRDS(here("data", "raw"#, "Survival"
                       , "survival-data.Rds"))

str(s.data)

survival_boxplot_summer <- ggplot(s.data %>% filter(Season == "Summer")
                                  , aes(x = Subspecies
                                        , y = Time
                                        , fill = Task)) +
  geom_boxplot(size = 0.7)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  #coord_flip()+
  facet_grid(cols = vars(Temperature)
             , rows = vars(Silica)
             #, scales = "free"
             )+
  theme_classic()+
  theme(strip.text = element_text(face = "bold"
                                  , size = 12)
        , strip.background = element_blank()
        , legend.key.size = unit(5, "mm")
        , legend.spacing.y = unit(0.2, "lines")
        , legend.title = element_text(size = 10)
        , legend.text = element_text(size = 9)
        , axis.text = element_text(size = 9)
        , axis.text.x = element_text(size = 9
                                     , face = "italic")
        , axis.title = element_text(size = 10))+ 
  labs(title = "Survival of summer bees"
       , y = "Survival time (min)")

survival_boxplot_winter1 <- ggplot(s.data %>% 
                                    filter(Season == "Winter")
                                  , aes(x = Subspecies
                                        , y = Time
                                        , fill = Task))+
  geom_boxplot(size = 0.7)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  #coord_flip()+
  facet_grid(cols = vars(Temperature)
             , rows = vars(Silica)
             #, scales = "free"
  )+
  theme_classic()+
  theme(strip.text = element_text(face = "bold"
                                  , size = 12)
        , strip.background = element_blank()
        , legend.key.size = unit(5, "mm")
        , legend.spacing.y = unit(0.2, "lines")
        , legend.title = element_text(size = 10)
        , legend.text = element_text(size = 9)
        , axis.text = element_text(size = 9)
        , axis.text.x = element_text(size = 9
                                     , face = "italic")
        , axis.title = element_text(size = 10))+ 
  labs(title = "Survival of winter bees"
       , y = "Survival time (min)")

survival_boxplot_winter2 <- ggplot(s.data %>% 
                                     filter(Season == "Winter" &
                                            !Collection.Hive == 3
                                     )
                                   , aes(x = Subspecies
                                         , y = Time
                                         , fill = Task))+
  geom_boxplot(size = 0.7)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  #coord_flip()+
  facet_grid(cols = vars(Temperature)
             , rows = vars(Silica)
             #, scales = "free"
  )+
  theme_classic()+
  theme(strip.text = element_text(face = "bold"
                                  , size = 12)
        , strip.background = element_blank()
        , legend.key.size = unit(5, "mm")
        , legend.spacing.y = unit(0.2, "lines")
        , legend.title = element_text(size = 10)
        , legend.text = element_text(size = 9)
        , axis.text = element_text(size = 9)
        , axis.text.x = element_text(size = 9
                                     , face = "italic")
        , axis.title = element_text(size = 10))+ 
  labs(title = "Survival of winter bees"
       , y = "Survival time (min)")
