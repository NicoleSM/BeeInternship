# This script contains the code to perform the analyses illustrated in the 
# Cluster Analysis section of the guide. The code still uses "Task"
# as the factor for the comparisons (Just as in the guide).
# DO NOT FORGET TO APROPIATELY CHANGE THIS AS NEEDED FOR YOUR OWN ANALYZES
# (PROBABLY WILL NOT WORK OTHERWHISE).

# Install pacman package if necessary ----
#install.packages("pacman")

# Load required packages ----
NPacks <- c("here","vegan","ggdendro", "tidyverse")
pacman::p_load(char=NPacks)

# Load the data frames ----
load(here("data", "processed", "Flying", "data-frames.Rdata"))


master.datenclust <- cbind.data.frame(group = as.factor(grouping_info$Task)
                                      , master.daten
                                      , row.names = paste(grouping_info$Task
                                                          , grouping_info$Individual))
str(master.datenclust)

ddata <- vegan::vegdist(master.datenclust %>% select(-group)
                        , method = "bray") %>% 
  hclust(method  = "ward.D2") %>% 
  as.dendrogram %>% 
  dendro_data()

dlabs <- label(ddata)
row.names(dlabs) <- dlabs$label
dlabs <-  merge(dlabs, master.datenclust[1], by = "row.names", all = T)

cluster_p <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_text(data = dlabs, aes(label = label, x = x, y = y-0.01, color = group
                              , hjust = "left"), size = 3.1) +
  scale_color_viridis_d()+
  coord_cartesian(ylim = c(-0.5, 1.5), xlim = c(0, 28))  +
  coord_flip(clip = "off") + 
  scale_y_reverse() + 
  theme(legend.position = "left"
        , legend.key = element_rect(fill = "white", color = "white")
        , panel.background = element_rect(fill = "white")
        , plot.margin = unit(c(1, 10, 1, 1), "lines")
        , axis.text.y = element_blank()
        , axis.ticks = element_blank()
        , axis.title = element_blank())
cluster_p
ggsave(here("figs", "Flying", "02_cluster-plot.png")
       , width = 30
       , height = 25
       , units = 'cm'
       , dpi = 'retina')

# End ----
## Report session information
capture.output(sessionInfo(), file = here("output", "Flying", "sInf_Script06.txt"))
# 
# ## Detach/unload packages
# lapply(NPacks, unloadNamespace)
# 
# ## Clear environment
# rm(list=ls())

