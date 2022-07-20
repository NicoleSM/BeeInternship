

# Instal the required packages if necessary ----

#install.packages(c('devtools','vegan'))
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis" 

                #If does not work, specify the library directory path
                 #, lib = "library directory path" 
 #            ) 

# Load the required packages ----
NPacks<-c("here", "vegan", "pairwiseAdonis", "tidyverse")

pacman::p_load(char=NPacks)

# Load the data frames ----
load(here("AMelliMelli-CHC-data", "data-frames.Rdata"))

## Defines the number of permutations (REQUIRED!!!!)
PERM<-how(nperm=999)

# Evaluate the difference in dispersion (average multivariate distance to median) between groups ----

##Distance based dispersion test
##Shows the average distance to median of each factor
BDTsk <- vegdist(master.daten, method="bray") %>% betadisper(grouping_info$Task)
BDTsk

## also Eigenvalue based method
boxplot(BDTsk, xlab="Task performance")

##F Test
##Permutation test
###If p-value is <0.05 difference in dispersion is significative
###If p-value is >0.05 dispersion is equal
set.seed(123456)
testBDTsk <- permutest(BDTsk, 
                       pairwise = F, # Set it as TRUE if the factor contains more than 2 groups
                       permutations = PERM)
testBDTsk
tBDTsk <- permustats(testBDTsk)
summary(tBDTsk) # SES and conf. interval for details check ?permutest()
densityplot(tBDTsk)
qqmath(tBDTsk)

# Perform the PERMANOVA
set.seed(123456)
PMVTsk <- adonis2(master.daten ~ Task, data = grouping_info, permutations = PERM, method = 'bray')
PMVTsk
tPMVTsk <- permustats(PMVTsk)
summary(tPMVTsk)
densityplot(tPMVTsk)
qqmath(tPMVTsk)

# Citation ----
citation("vegan")
RStudio.Version()
citation() #To cite R