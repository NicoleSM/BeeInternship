# This script contains the code to perform the analyses illustrated in the 
# Non-metric Multidimensional Scaling section of the guide. 
# The code still uses "Task" as the factor for the comparisons (Just as in the guide).
# DO NOT FORGET TO APROPIATELY CHANGE THIS AS NEEDED FOR YOUR OWN ANALYZES
# (PROBABLY WILL NOT WORK OTHERWHISE).

# Install the necessary packages if the are not yet ----
NPacks<-c("here", 'vegan', 'goeveg', "plotly", "viridis", "tidyverse")
#install.packages(NPacks)

# Load required packages ----
pacman::p_load(char=NPacks)

# Load the data frames ----
load(here("AMelliMelli-CHC-data", "processed", "data-frames.Rdata"))


# Perform the NMDS ----
## Determine the number of dimensions for the NMDS
### As a rule of thumb, an NMDS ordination with a stress value around or above 0.2 is deemed suspect 
### and a stress value approaching 0.3 indicates that the ordination is arbitrary. Stress values equal
### to or below 0.1 are considered fair, while values equal to or below 0.05 indicate good fit
dimcheckMDS(master.daten, distance = "bray", trymax=400)

## Builds the NMDS
set.seed(123456) #Always set the seed before running the metaMDS comand
sol<-metaMDS(master.daten, distance = "bray",
              k = 2, #For 2D NMDS
              try = 200,  
              #k = 3, #For 3D NMDS
             trymax = 400)
sol

## Evaluating the NMDS
stressplot(sol)

# Makes a new data frame, which will be used for plotting
# note: Be sure to define the number of dimensions properly or it will not work
NMDS <- data.frame(NMDS1 = sol$point[, 1], NMDS2 = sol$point[, 2]
                   #, NMDS3 = sol$point[, 3]  # Only if the NMDS is with 3D
                   , grouping_info)
NMDS

# 2D NMDS ----

# Centroids
centTsk <- aggregate(cbind(NMDS1, NMDS2) ~ Task, data = NMDS, FUN = mean)
centTsk

# Ordispider segments 
spiderTsk<-merge(cbind(NMDS[1:2], 'Task' = NMDS$Task)
                 , setNames(centTsk,c('Task', 'oNMDS1', 'oNMDS2'))
                 , by = 'Task', sort = F)
spiderTsk

nmdsPlot<-ggplot(data = NMDS,aes(x = NMDS1, y = NMDS2, color = Task)) +
  stat_ellipse(#aes(linetype = Task), 
    size = 0.72, level = 0.95) +
  geom_segment(data = spiderTsk
               , aes(xend = oNMDS1, yend = oNMDS2#, linetype = Task
               )#, alpha = 0.75
               , size = 0.72) + 
  geom_point(aes(shape = Task, fill = Task), alpha = 0.75, size = 3) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c("In-hive workers" = 25
                                , "Out-hive workers" = 23)) +
  #scale_linetype_discrete() +
  coord_equal() +
  theme_classic()
nmdsPlot
ggsave(here("figs", "03_nmds-plot.png")
       , width = 30
       , height = 25
       , units = 'cm'
       , dpi = 'retina')

# Plot 3D NMDS ----

# ## Establish the number of groups to be compared in the plot
# ## It is required for the color scale
# num.groups <- length(unique(NMDS.3d$Task))
# 
# ## Define the scene's set-up, establishing the position of the camera
# scene = list(camera = list(eye = list(x = -2, y = 1.35, z = 0.5))) 
# 
# nmdsPlot <- plot_ly(NMDS, x = ~NMDS1, y = ~NMDS3, z = ~NMDS2
#                     , color = ~Task
#                     , colors = viridis_pal(option = "D")(num.groups)
#                     , alpha = 1, size = 20, sizes = c(100, 100)
#                     , symbol = ~Task
#                     , symbols = c("x", "diamond"))  %>%
#   layout(scene = scene)
# nmdsPlot

# End ----
## Report session information
capture.output(sessionInfo(), file = here("output", "sInf_Script07.txt"))
# 
# ## Detach/unload packages
# lapply(NPacks, unloadNamespace)
# 
# ## Clear environment
# rm(list=ls())

