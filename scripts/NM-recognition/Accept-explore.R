
# Install and/or load packages ----
NPacks <- c("skimr","GGally", "ggplot2", "tidyr", "dplyr", "here")

#install.packages(c(NPacks,"pacman"))

## Load all the packages and try to install them if they are not (Uses/require pacman package)
pacman::p_load(char = NPacks)

# Load the acceptance data ----
accept_daten <- readRDS(here("data", "processed"
                             , "NM-Recognition", "accept_daten.Rds"))
str(accept_daten)

#----
accept_daten %>% group_by(location, Subspecies) %>% skim()

png(here("figs", "NMRec_02-nice.pairplot.png")
    , width = 800, height = 800, units = "px")
ggpairs(accept_daten %>% select(!n_nestmates & !n_non.nestmates)
        , aes(color = Subspecies, alpha = 0.8)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d()+
  theme_classic()
dev.off()

ggplot(accept_daten
       , aes(x = Subspecies
             , y = nestmates
             , fill = Subspecies)) +
 # ggdist::stat_slab(point_color = NA
  #                  , trim = F
   #                 , justification = -0.4) +
  geom_boxplot(aes(x = Subspecies)
               , width = 0.3
               , outlier.shape = NA) +
  geom_point(shape = 21
             , size = 2.3
             , alpha = 0.8
             , position = position_jitter(width = 0.1
                                          , height = 0
                                          , seed = 1)) +
  coord_flip() +
  scale_fill_viridis_d() +
  facet_wrap(vars(location)) +
  theme_classic()


# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output"
                             , "sInf_NM_Recognition-script_02.txt"))

## Detach/unload packages
lapply(NPacks, unloadNamespace)