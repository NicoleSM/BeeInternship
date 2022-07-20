
# Install and/or load packages ----
NPacks <- c("skimr","GGally", "ggplot2", "tidyr", "dplyr", "here")

#install.packages(c(NPacks,"pacman"))

## Load all the packages and try to install them if they are not (Uses/require pacman package)
pacman::p_load(char = NPacks)

# Load the acceptance data ----
nm.rec_daten <- readRDS(here("data", "processed"
                             , "NM-Recognition", "nm.rec_daten.Rds"))
str(nm.rec_daten)

factrs_omit <- c("Bee.number"
                 , "Fostered"
                 , "Foster.hive"
                 , "group"
                 , "Acceptance")
# ggpairs(nm.rec_daten,
#         columns = nm.rec_daten %>% select(-all_of(factrs_omit)) %>% colnames()
#         , aes(color = group)
#         , legend = c(1, 1)
#         , cardinality_threshold = 16) +
#   scale_color_viridis_d() +
#   scale_fill_viridis_d() +
#   facet_wrap(vars(location)) +
#   theme_classic()
# ggsave("")

nm.rec_hist <- ggplot(nm.rec_daten) +
  geom_histogram(aes(x = Acceptance, fill = group)
                 , stat = "count"
                 , position = position_dodge()
                 #, alpha = 1
                 ) +
  scale_x_discrete(labels = c("Rejected", "Accepted")) +
  scale_fill_viridis_d(option = "inferno"
                       , labels = c("Nestmate", "Non-nestmate")) +
  facet_grid(cols = vars(Subspecies)
             , rows = vars (location)
             , scales = "free_y") +
  theme_classic() +
  theme(legend.position = "left"
        , strip.background = element_blank()
        , strip.text.x = element_text(face = "italic")
        , strip.text = element_text(face = "bold"
                                    , size = 10)
        , legend.title = element_text(size = 10)
        , legend.text = element_text(size = 8.7)
        , axis.title.x = element_blank()
        , axis.text = element_text(size = 8.7))

# ggplot(nm.rec_daten, aes(x = group
#                          , y = Acceptance
#                          , fill = Subspecies)) +
#   geom_point(shape = 21
#              , size =  4
#              , alpha = 0.8
#              , position = position_jitter(width = 0.2
#                                           , height = 0.08
#                                           , seed = 1)) +
#   scale_fill_viridis_d() +
#   facet_grid(cols = vars(Subspecies)
#              , rows = vars (location)
#              , scales = "free_x") +
#   theme_classic()

