# Install the necessary packages ----

  #Vector holding the list of packages that will be used
NPacks<-c("here", "reshape2", "skimr","MVN", "ggdist", "GGally", "ggpubr"
          , "dunn.test", "vegan", "tidyverse")

#install.packages(NPacks)

# Load required packages ----
pacman::p_load(char=NPacks)

# FUNCTIONS ----
## In this section the functions that will be used latter in the script are defined
## Without running the code on this section in advance, some parts of the script
## will not work as they rely on these functions

## Data frame by hydrocarbon classes ----
cclasses_df <- function(daten, grouping.info, Comps.data, fuse.methyls=F){
  require(dplyr)
 
  # Data set per hydrocarbon class
  ## Vector to store compund classes` names
  cclasses <- c()
  
  for (c.class in unique(Comps.data$Class)) {
    cclasses <- c(cclasses, c.class) 
    tmp <- as.data.frame(t(daten)) %>% 
      filter(Comps.data$Class == c.class)
    assign(c.class, tmp)
  }
  
  Sampls_index <- data.frame('Sample name' = rownames(daten)
                             , rank = rank(grouping.info$Individual)
                             , stringsAsFactors  =  FALSE)
  
  if(fuse.methyls == T){
    Prop.CompsClass <- data.frame(grouping.info
                                  , Alkanes = colSums(Alkane)
                                  , Alkenes = colSums(Alkene)
                                  , Alkadienes = colSums(Alkadiene)
                                  , Methyls = rbind(Methyl
                                                    , Dimethyl) %>%
                                    colSums()
                                  , row.names  =  Sampls_index$Sample.name)
  } else{
    # Has to be adjusted to let the function work
    Prop.CompsClass <- data.frame(grouping.info
                                  , Alkanes = colSums(Alkane)
                                  , Alkenes = colSums(Alkene)
                                  , Alkadienes = colSums(Alkadiene)
                                  , Monomethyls = colSums(Methyl)
                                  #, Dimethyls = colSums(Dimethyl)
                                  , row.names  =  Sampls_index$Sample.name)
  }
  
  Prop.CompsClass %>% as_tibble()
}

## Data frame of peaks richness per hydrocarbon class ----
richness_cclasses <- function(daten, grouping.info, Comps.data, faktor){
  Alka <- as.data.frame(t(daten)) %>%
    filter(Comps.data$Class=='Alkane')
  Alke <- as.data.frame(t(daten)) %>%
    filter(Comps.data$Class=='Alkene')
  Alkadi <- as.data.frame(t(daten)) %>%
    filter(Comps.data$Class=='Alkadiene')
  Met <- as.data.frame(t(daten)) %>%
    filter(Comps.data$Class=='Methyl')
  Dime <- as.data.frame(t(daten)) %>%
    filter(Comps.data$Class=='Dimethyl')
  # Trime <- as.data.frame(t(daten)) %>%
  #   filter(Comps.data$Class=='Trimethyl')
  # Tetrame <- as.data.frame(t(daten)) %>%
  #   filter(Comps.data$Class=='Tetramethyl')
  
  AlkaRows <- c()
  Alkanes <- data.frame()
  AlkeRows <- c()
  Alkenes <- data.frame()
  AlkadRows <- c()
  Alkadienes <- data.frame()
  MetRows <- c()
  Methyls <- data.frame()
  DimRows <- c()
  Dimethyls <- data.frame()
  
  for(n in levels(grouping.info %>% pull(faktor))){
    tAlka <- as.data.frame(t(Alka)
                           , row.names = rownames(as.data.frame(t(Alka))))
    tAlka <- tAlka %>%
      filter(as.factor(grouping.info[, faktor]) == n) %>% 
      colSums()
    tAlka[tAlka[] > 0] <- 1
    tAlka <- sum(tAlka)
    Alkanes <- rbind(Alkanes, tAlka = sum(tAlka))
    AlkaRows <- c(AlkaRows, n)
    
    tAlke <- as.data.frame(t(Alke)
                           , row.names = rownames(as.data.frame(t(Alke))))
    tAlke <- tAlke %>%
      filter(as.factor(grouping.info[, faktor]) == n) %>% 
      colSums()
    tAlke[tAlke[] > 0] <- 1
    tAlke <- sum(tAlke)
    Alkenes <- rbind(Alkenes, tAlke = sum(tAlke))
    AlkeRows <- c(AlkeRows, n)
    
    tAlkad <- as.data.frame(t(Alkadi)
                            , row.names = rownames(as.data.frame(t(Alkadi))))
    tAlkad <- tAlkad %>%
      filter(as.factor(grouping.info[, faktor]) == n) %>%
      colSums()
    tAlkad[tAlkad[] > 0] <- 1
    tAlkad <- sum(tAlkad)
    Alkadienes <- rbind(Alkadienes
                        , tAlkad = sum(tAlkad))
    AlkadRows <- c(AlkadRows, n)
    
    tMet <- as.data.frame(t(Met)
                          , row.names = rownames(as.data.frame(t(Met))))
    tMet <- tMet %>% 
      filter(as.factor(grouping.info[, faktor]) == n) %>% 
      colSums()
    tMet[tMet[] > 0] <- 1
    tMet <- sum(tMet)
    Methyls <- rbind(Methyls
                     , tMet = sum(tMet))
    MetRows <- c(MetRows, n)
    
    tDim <- as.data.frame(t(Dime)
                          , row.names = rownames(as.data.frame(t(Dime))))
    tDim <- tDim %>%
      filter(as.factor(grouping.info[, faktor]) == n) %>% 
      colSums()
    tDim[tDim[] > 0] <- 1
    tDim <- sum(tDim)
    Dimethyls <- rbind(Dimethyls, tDim = sum(tDim))
    DimRows <- c(DimRows, n)
  }
  
  rownames(Alkanes) <- AlkaRows
  colnames(Alkanes) <- "Alkanes"
  rownames(Alkenes) <- AlkeRows
  colnames(Alkenes) <- "Alkenes"
  rownames(Alkadienes) <- AlkadRows
  colnames(Alkadienes) <- "Alkadienes"
  rownames(Methyls) <- MetRows
  colnames(Methyls) <- "Methyls"
  rownames(Dimethyls) <- DimRows
  colnames(Dimethyls) <- "Dimethyls"
  
  tmp <- data.frame(group = as.factor(row.names(Alkanes))
                    , Alkanes
                    , Alkenes
                    , Alkadienes
                    , Methyls
                    , Dimethyls)
  tmp <- tidyr::pivot_longer(tmp, cols = !where(is.factor)
                             , names_to = "variable"
                             , values_to = "value") %>% 
    as.data.frame()
  tmp$variable <- factor(tmp$variable, levels = unique(tmp$variable))
  
  CClasses <- levels(tmp$variable)
  CClasses[][CClasses[] == "Methyls"] <- "Monomethyl alkanes"
  CClasses[][CClasses[] == "Dimethyls"] <- "Dimethyl alkanes"
  
  levels(tmp$variable) <- CClasses
  
  tmp
}

## Chain length data set ----
## Function to generate a data set with the weighted average chain length per
## individual
cl_df <- function(daten, Comps.data, grouping.info){
  cl.df <- data.frame(row.names = rownames(daten))
  Ncols <- as.character()
  for (i in unique(Comps.data$Chain.length)) {
    #We need to sum the abundance of the compounds#
    Tmp <- as.data.frame(t(daten)) %>%
      filter(Comps.data$Chain.length == i) %>%
      colSums() %>% 
      as.numeric()
    cl.df <- cbind.data.frame(cl.df, Tmp)
    Ncols <- c(Ncols, i)
  }
  rm(i, Tmp)
  colnames(cl.df) <- Ncols
  rm(Ncols)
  
  cl.df <- cbind(grouping.info, cl.df)
  
  merge(grouping.info
        , cl.df %>% 
          pivot_longer(cols = !where(is.factor)
                       , names_to = "Chain.length"
                       , values_to = "Abundance") %>%  
          mutate(Chain.length = as.numeric(Chain.length)) %>%
          group_by(Individual) %>%
          summarise(cl_w.mean = weighted.mean(Chain.length
                                              , Abundance))
        , by = "Individual") %>% as_tibble()
}

# Load the data frames ----
load(here("AMelliMelli-CHC-data", "processed", "data-frames.Rdata"))

## Creates a data frame with the abundance of each compounds class per sample and the factors to group the samples
Prop_CompsClass <- cclasses_df(master.daten
                               , grouping_info
                               , master.Comps)

head(Prop_CompsClass)
str(Prop_CompsClass)


# Peaks per hydrocarbon class  ----
# Table with the count of compounds per class of each group
rich_task <- richness_cclasses(master.daten
                  , grouping_info
                  , master.Comps
                  , faktor = "Task")
rich_task

ggplot(rich_task, aes(x = variable, y = value, fill = group)) +
  geom_bar(position = position_dodge(), stat = "identity", alpha = 1) +
  scale_fill_viridis_d() +
  xlab("Class of hydrocarbons") + ylab("No. of peaks") + labs(fill = "Task") +
  theme_classic()

# Hydrocarbon classes' abundance per group ----

# summary stats and plots

## Pair plots summarizing the data set
pdf('.pdf'
    ,width = 28, height = 28)
PP<-Prop_CompsClass[-4]
head(PP)
ggpairs(PP, cardinality_threshold = NULL, aes(color=Task, alpha=0.7)) + 
  scale_color_viridis_d() +
  scale_fill_viridis_d() + 
  theme_classic()
dev.off() # Necessary every time you plot after using the bmp() or pdf() function

## Detailed summary
Prop_CompsClass %>% group_by(Task) %>% skim()

## Normality test (NOT complete YET, but working)----
##NOTE: GRaphic methods may fail showing an error
# Most probably due to colineality, in that case you can be pretty sure the data tested are
## normal on the level (multivar/univar/group/Compound class) that you are testing for it
### Multivariate normality (Most probably will not occur, but let's follow the formality)

#### A data frame grouping the data for each factor
CClass.Task<-as.data.frame(Prop_CompsClass[5:ncol(Prop_CompsClass)]) %>% 
  mutate(group=as.factor(Prop_CompsClass$Task)) %>%
  as.data.frame()
str(CClass.Task)

### Now the multivariate-normality test
### If you get: Error in solve.default(S, tol = tol) : Lapack routine dgesv: system is exactly singular: U[#,#] = 0
### Most probably there are colinear variables in at least one group, especially if one or more compounds are absent in a group
### At this point you can be very confident that there is no multivariate normality

pdf('CompsClassNFQQplots.pdf',width = 5, height = 5)
mvn(CClass.Task, subset = "group"
    , mvnTest = "hz"
    , multivariatePlot = "qq"
    , desc = FALSE
    , univariateTest="Lillie"
    , univariatePlot = "histogram")
dev.off() # Necessary every time you plot after using the bmp() or pdf() function


## Homocedasticity test----

CClass.Task %>% 
  select(!where(is.factor)) %>% 
  lapply(bartlett.test, CClass.Task$group)

## Statistical tests ----
CClasses <- Prop_CompsClass %>% 
  select(!where(is.factor)) %>% 
  colnames()

### Parametric ----
#### t-test----
#### Assumes normality and variance homocedasticity. Alternative: Wilcoxon rank-sum/Mann-Whitney U
#### Cannot deal with factors of more than two-levels. Alternative: ANOVA
for (i in CClasses) {
  print(paste("i", i, sep = " = "))
  test <- t.test(data=CClass.Task, CClass.Task[,i] ~ group)
  print(test)
  assign(paste("test", i, sep = "_"), test)
  rm(test, i)
}

#### ANOVA----
#### Assumes normality and Homocedasticity of variances. Alternative: Kruskal-Wallis
for (i in CClasses) {
  print(paste("i", i, sep = " = "))
  aov <- aov(CClass.Task[, i] ~ group, CClass.Task)
  print(aov)
  sum_aov <- summary(aov)
  print(sum_aov)
  ph <- TukeyHSD(aov) # Post-hoc Tukey test
  print(ph)
  test <- list(aov, sum_aov, ph)
  assign(paste("test", i, sep = "_"), test)
  rm(test, i, aov, sum_aov, ph)
}

### Non-parametric ----
#### Wilcoxon rank-sum/Mann-Whitney U----
for (i in CClasses) {
  print(paste("i", i, sep = " = "))
  test <- wilcox.test(data = CClass.Task, CClass.Task[, i] ~ group)
  print(test)
  assign(paste("test", i, sep = "_"), test)
  rm(test, i)
}

#### Kruskal-Wallis rank-sum----
for (i in CClasses) {
  print(paste("i", i, sep = " = "))
  kw <- kruskal.test(data = CClass.Task, CClass.Task[, i] ~ group)
  print(kw)
  ph <- dunn.test(CClass.Task[, i]  # Post-hoc Dunn's test
                             , CClass.Task$group
                             , kw=FALSE
                             , table=FALSE
                             , list=TRUE
                             , method = 'bonferroni')
  print(ph)
  test <- list(kw, ph)
  assign(paste("test", i, sep = "_"), test)
  rm(test, i, kw, ph)
}

## Plots ----

# Violin plots
for(i in CClasses){
  YLAB<-paste(i,"(%)",sep = " ")
  p<-ggplot(Prop_CompsClass
            , aes(x = Task, y = Prop_CompsClass %>% pull(i), fill = Task)) +
    geom_violin(color = "white", alpha=0.7, trim = F) +
    geom_boxplot(width = 0.2, color = "grey10", outlier.shape = NA) +
    geom_point(aes(fill=Task), color = "grey30", alpha = 0.75
               , shape = 21, stroke = 1.1, size = 1.3
               , position = position_jitter(width = 0.08, seed = 1)) +
    scale_fill_viridis_d() +
    theme_classic() +
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
    labs(y=YLAB)
  print(p)
  rm(i,YLAB,p)
}

# Rain cloud / Little prince plots
for(i in CClasses){
  YLAB<-paste(i
              ,"(%)"
              ,sep = " ")
  p<-ggplot(Prop_CompsClass
            , aes(x = Task
                  , y = Prop_CompsClass %>% 
                    pull(i)
                  , fill = Task)) +
    ggdist::stat_slab(width = 0.4
                      , alpha = 0.7
                      , point_color = NA
                      # , slab_type = "histogram"
                      # , breaks = 5
                      , trim = F
                      , justification = -0.7) +
    geom_boxplot(aes(x=as.numeric(Task)+0.25)
                 , width = 0.2
                 , color = "grey10"
                 , outlier.shape = NA) +
    geom_point(aes(fill=Task)
               , color = "grey30"
               , alpha = 0.75
               , shape = 21
               , stroke = 1.1
               , size = 1.3
               , position = position_jitter(width = 0.08
                                            , height = 0
                                            , seed = 1)) +
    scale_fill_viridis_d() +
    coord_flip() +
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(y=YLAB)
  print(p)
  rm(i,YLAB,p)
}

# Chain-length analysis ----
## Relative abundance of equivalent chain length ####
Prop_chain.length <- cl_df(master.daten, master.Comps, grouping_info)

# Little prince plot
Prop_chain.length %>% 
  ggplot(aes(x = Task
             , y = cl_w.mean
             , fill = Task)) +
  ggdist::stat_slab(width = 0.73
                    , alpha = 0.7
                    , point_color = NA
                    , trim = F
                    , justification = -0.32) +
  geom_boxplot(aes(x=as.numeric(Task)+0.21)
               , width = 0.2
               , color = "grey10"
               , outlier.shape = NA) +
  geom_point(color = "grey30"
             , alpha = 0.75
             , shape = 21
             , stroke = 1.1
             , size = 1.3
             , position = position_jitter(width = 0.08
                                          , height = 0
                                          , seed = 1)) +
  # geom_boxplot(color = "grey10"
  #              , alpha = 0.8) +
  scale_fill_viridis_d()+ 
  coord_flip() +
  theme_classic() +  
  theme(legend.key.size = unit(5, "mm")
        , legend.spacing.y = unit(0.2, "lines")
        , legend.title = element_text(size = 10)
        , legend.text = element_text(size = 8)) +
  labs(y = "Weighted mean of hydrocarbons chain length")

# Citation of used packages----
lapply(NPacks, citation)
citation() #To cite R
