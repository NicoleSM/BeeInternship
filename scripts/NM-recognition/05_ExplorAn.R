# This script contains the code to perform the analyses illustrated in the 
# Exploratory Data Analysis section of the guide. The code still uses "Task"
# as the factor for the comparisons (Just as in the guide).
# DO NOT FORGET TO APROPIATELY CHANGE THIS AS NEEDED FOR YOUR OWN ANALYZES
# (PROBABLY WILL NOT WORK OTHERWHISE).

# Install pacman package if necessary ----
#install.packages("pacman")

# Load required packages ----
NPacks<-c("here", "GGally", "tidyverse", "skimr", "MVN")
pacman::p_load(char=NPacks)

# Load the data frames ----
load(here("AMelliMelli-CHC-data", "processed", "data-frames.Rdata"))

# Data frame with factors information
DEA<- cbind.data.frame(grouping_info, master.daten)
str(DEA)

# Visualize the data set ----
## Pair plot ####
##  If your data set is too big, this might not be a good option

### Base plot version

### Creates and export to the project folder a PDF with a plot that allows
### initial visual exploration of the data set structure
# pdf(here("figs", "nice.pairplot.pdf"),width = 150, height = 150)
# plot(DEA %>% select(-Individual)) #The column with the individual number is excluded
# dev.off() # Run this as soon as the plot process finish to store the pdf

### GGAlly version
### Play with the size of the PDF to assure good visualization
### CHC data sets tend to be too big, so the utility of this plot
### is very limited the resulting plot can be too overwhelming

# Creates and export ta PDF with the nicer pair plot
# pdf(here("figs", "nice.pairplot.pdf"), width = 150, height = 150)  
# ggpairs(DEA %>% select(-Individual) # Exclude column with individual number  
#         , cardinality_threshold = NULL  
#         # Colors to differ among groups of a selected factor (e.g. Task)
#         ,aes(color = Task 
#              ,alpha = 0.6)) +
#   scale_color_viridis_d() +
#   scale_fill_viridis_d()+
#   theme_classic()
# dev.off()  # Run this as soon as the plot process finish to store the pdf

## Heatmap ####
melt_DEA <- pivot_longer(DEA, cols = !where(is.factor)
                                , names_to = "Compound"
                                , values_to = "Abundance")
melt_DEA$Compound <- factor(melt_DEA$Compound
                            , levels = unique(melt_DEA$Compound))
melt_DEA <- melt_DEA %>% as.data.frame()
str(melt_DEA)

ggplot(melt_DEA, aes(x = Individual, y = Compound)) +
  geom_raster(aes(fill = log1p(Abundance))) +
  scale_fill_viridis_c(option = "turbo") +
  # Y axis is sorted from shorter chain length to longer chain length
  # labels correspond to the compound name
  scale_y_discrete(limits = rev
                   , labels = master.Comps$Compound %>%
                     rev()) +
  # Produce facets wrapping the plot regarding the task groups
  facet_wrap(vars(Task), scales = "free_x") + 
  theme_classic()
ggsave(here("figs", "01_heatmap.png")
       , width = 30
       , height = 25
       , units = 'cm'
       , dpi = 'retina')

# Summarize the the data set ----
DEA %>% group_by(Task) %>% skim() %>% print()

# Test for normality ----

# A data frame grouping the data with regards to the selected factor
normal.task <- as.data.frame(DEA[5:ncol(DEA)]) %>%
  mutate(group = as.factor(DEA$Task)) %>%
  as.data.frame()

# Now the normality test
# If you get: 
# Error in solve.default(S, tol = tol) : Lapack routine dgesv: system is exactly singular: U[#,#] = 0
# This means that there are colinear variables in at least one group, 
# especially if one or more compounds are absent in a group.
# At this point you can be very confident that there is no multivariate normality
# mvn.task<-mvn(normal.task, subset = "group" 
#                    , mvnTest = "hz"
#                    , desc=TRUE
#                    , univariateTest="Lillie")
# mvn.task

# End ----
## Report session information
capture.output(sessionInfo(), file = here("output", "sInf_Script05.txt"))
# 
# ## Detach/unload packages
# lapply(NPacks, unloadNamespace)
# 
# ## Clear environment
# rm(list=ls())
