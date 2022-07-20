# Install the necessary packages ----
NPacks<-c("pacman", "here" ,"randomForest", "dplyr", "ggplot2")
#install.packages(NPacks)

# Load required packages ----
pacman::p_load(char=NPacks)

# Load the data frames ----
load(here("AMelliMelli-CHC-data", "data-frames.Rdata"))

tsk.daten=data.frame(treatment=as.factor(grouping_info$Task), master.daten)
tsk.daten

# Creates the randomForest ----
## If classification between groups of a factor will be compared among groups of another factor,
##  thus using subsets of the data set, then a randomforest classification with the same parameters
##  (seed, mtry, importance, do.trace, ntree, proximity, keep.forest, etc.) must be performed for
##  each subset's data frame.

length(colnames(master.daten))
set.seed(92) # Set it everytime before building the random forest to be able to replicate results
tsk.rf <- randomForest(treatment~., tsk.daten,
                       mtry = sqrt(length(colnames(master.daten)))
                       ,importance=TRUE,
                       do.trace=1000, ntree=10000, proximity=TRUE,
                       keep.forest=TRUE)
tsk.rf
plot(tsk.rf)
#note: rule for mtry: mtry=sqrt(n_compounds)plot(tsk.rf) 

# MDS plot with samples: Shows the relationship between them ----
## This is not an NMDS (This one is metric) it only shows the relationships between groups found by the Random Forest classification
## Status refers to the treatment: Remember to establish the same treatment that was used to build the Random Forest
## or the resulting plot will be ABSOLUTELY wrong
mds.stuff <- dist(1-tsk.rf$proximity) %>% cmdscale(eig = TRUE, x.ret = TRUE)
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample = rownames(mds.values)
                       , X = mds.values[, 1]
                       , Y = mds.values[, 2]
                       , Indv = grouping_info$Individual
                       , Status=tsk.daten$treatment
)

png(filename = 'SvWRF-MDS.png',width = 800, height = 800, units = "px", res = NA)
#pdf('RF-MDS.pdf',width = 15, height = 15)
ggplot(data=mds.data, aes(x=X, y=Y, label= Indv)) +
  geom_text(aes(x=X, y=Y, color=Status), cex = 4, fontface = "bold") +
  scale_color_viridis_d() +
  theme_minimal() +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep = "")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep = "")) +
  #coord_equal() +
  ggtitle("MDS plot using (1 - Random Forest proximities)")
dev.off()

# Variable Importance ----
## Store the importance of the variables in a data frame
## The data frame holds the Mean Decrease of Accuracy for each group and for the overall classification
var_importance <- randomForest::importance(tsk.rf, scale = TRUE) %>% as.data.frame()
cnames <- colnames(master.daten)
rownames(var_importance) <- cnames
var_importance <- var_importance %>% arrange(desc(var_importance$MeanDecreaseAccuracy))
var_importance

## Exports the var_importance table as a CSV file, add the name of the factor(s) to the name of the file to identify it
write.csv(var_importance, file=here("AMelliMelli-CHC-data", "var_importanceSvW.csv"), sep=";", row.names = TRUE, col.names = TRUE)

## Variable importance plot
png(filename = 'VarImpPlotRF-SvW.png',width = 400, height = 600, units = "px", res = NA)
p <- ggplot(var_importance
            , aes(y = reorder(row.names(var_importance), MeanDecreaseAccuracy)
                  , x = MeanDecreaseAccuracy))
p <- p + ylab("Compound") +
  xlab("Variable Importance (Mean Decrease of Accuracy Index)")

p <- p + 
  geom_point(color = "#440154FF", size = 2) +
  geom_vline(xintercept = mean(var_importance$MeanDecreaseAccuracy)
             , linetype = "solid"
             , color = "grey50") +
  geom_vline(xintercept = quantile(var_importance$MeanDecreaseAccuracy
                                   , probs = c(0.25, 0.5, 0.75))
             , linetype = "dashed"
             , color = "grey50") +
  scale_colour_viridis_d()

p + theme_classic() +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.y = element_line(size = 0.5))
dev.off()

#Partial Dependence Plots (Explore more!!)
p<-partialPlot(tsk.rf,pred.data = WBseasn,x.var = C31,which.class = "Summer")
p1
p2<-partial(tsk.rf, pred.var="C31", which.class = "Summer"
            , type="auto", plot=T, rug=T, plot.engine ="ggplot2")
p2



# End ----
## Report session information
capture.output(sessionInfo(), file = here("output", "sInf_Data-frames-script_08.txt"))

## Detach/unload packages
lapply(NPacks, unloadNamespace)