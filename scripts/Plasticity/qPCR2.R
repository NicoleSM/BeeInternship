library(EasyqpcR)
library(reshape2)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(car)
library(FSA)
library(mdthemes)
library(ggpol)
library(matrixStats)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read data ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setwd("C:\\Users\\Nicole\\Documents\\BeeInternship\\data\\raw\\Plasticity\\data-analysis")
Standard <- read.csv("SK.csv", header=T, check.names = FALSE)
run1_2_CT <- read.csv("run1_2.csv",header=T,check.names = FALSE)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# calculate relative expression ----
## calculate efficiencies ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ef <- unlist(slope(data=Standard, q=c(50000000,
                                      5000000,
                                      500000,
                                      50000,
                                      5000,
                                      500),r=3, na.rm=TRUE))
Ef

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## calibrator factors ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CF1 <- calData(nrmData(data = run1_2_CT, #run1
                       r=3,
                       E=Ef,
                       Eerror=c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02),
                       nSpl=42,
                       nbRef=2,
                       Refposcol=1:2,
                       nCTL=5,
                       CF=c(1, 1, 1, 1, 1, 1),
                       CalPos=c(41:42),
                       trace = FALSE,
                       geo = TRUE,
                       na.rm = TRUE)[[3]][1:1,])

CF2 <- calData(nrmData(data = run1_2_CT, #run2
                       r=3,
                       E=Ef,
                       Eerror=c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02),
                       nSpl=42,
                       nbRef=2,
                       Refposcol=1:2,
                       nCTL=5,
                       CF=c(1, 1, 1, 1, 1, 1),
                       CalPos=c(41:42),
                       trace = FALSE,
                       geo = TRUE,
                       na.rm = TRUE)[[3]][2:2,])

############################################################################
# calculate Normalized RQ >>> Calibrated NRQ (for each run)
############################################################################
CNRQ1 <- nrmData(data = run1_2_CT,
                 r=3,
                 E=Ef,
                 Eerror=c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02),
                 nSpl=42,
                 nbRef=2,
                 Refposcol=1:2,
                 nCTL=5,
                 CF=CF1,
                 CalPos=c(41),
                 trace = FALSE,
                 geo = TRUE,
                 na.rm = TRUE)[[2]][c(1:20,41),]
CNRQ2 <- nrmData(data = run1_2_CT,
                 r=3,
                 E=Ef,
                 Eerror=c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02),
                 nSpl=42,
                 nbRef=2,
                 Refposcol=1:2,
                 nCTL=5,
                 CF=CF2,
                 CalPos=c(42),
                 trace = FALSE,
                 geo = TRUE,
                 na.rm = TRUE)[[2]][c(21:40,42),]

############################################################################
# aggregate all CNRQs
############################################################################
CNRQs <- rbind(CNRQ1,CNRQ2)
CNRQs

############################################################################
# auf die Kontrollgruppe normieren (geometrisches Mittel)
### Kontrollgruppe: Sommerbienen + Nicht-Heizer
############################################################################
ctlgroup <- as.matrix(CNRQs[c(1:5,22:26),])
ctlgeom <- colProds(ctlgroup)^(1/dim(ctlgroup)[1])
ctlgeom1 <- (as.data.frame(ctlgeom)[rep(1:(ncol(CNRQs)), each = nrow(CNRQs)), ])
ctlgeom2 <- as.data.frame(matrix(ctlgeom1, ncol = ncol(CNRQs), byrow = FALSE))
CNRQs_scaled_to_group <- CNRQs/ctlgeom2
CNRQs_scaled_to_group

############################################################################
# den Kalibrator entfernen
############################################################################
CNRQs_scaled_to_group <- CNRQs_scaled_to_group[-c(21,42), ]
CNRQs_scaled_to_group

CNRQs_scaled_to_group$subspecies <- as.factor(c(rep("carnica", 10),
                                        rep("iberiensis", 10),
                                        rep("carnica", 10),
                                        rep("iberiensis", 10)
))
CNRQs_scaled_to_group$role <- as.factor(c(rep("forager bee", 5),
                                  rep("nurse bee", 5),
                                  rep("forager bee", 5),
                                  rep("nurse bee", 5),
                                  rep("forager bee", 5),
                                  rep("nurse bee", 5),
                                  rep("forager bee", 5),
                                  rep("nurse bee", 5)
))

CNRQs_scaled_to_group
CNRQs_scaled_to_group_m <- melt(subset(CNRQs_scaled_to_group, select=c("subspecies","role","Des1","Des2","El1","El2")), id=c("subspecies","role"), variable.name = "gene", value.name = "expression")
CNRQs_scaled_to_group_m

ggplot(data = CNRQs_scaled_to_group_m, aes(x=subspecies,
                                   y=expression,
                                   fill=role
                                   )) +
  geom_boxjitter(jitter.shape = 21,
                 errorbar.draw = FALSE,
                 show.legend = TRUE) +
  facet_wrap(~gene,
             scales = "free"
             )+
  theme_classic()+
  labs(y='relative expression')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

