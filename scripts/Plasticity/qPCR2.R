library(EasyqpcR)
library(reshape2)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(magick)
library(cowplot)
library(ggsignif)
library(car)
library(FSA)
library(png)
library(mdthemes)
library(ggpol)
library(emmeans)
library(DHARMa)

 # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read data ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setwd("C:\\Users\\Nicole\\Documents\\BeeInternship\\data\\raw\\Plasticity\\data-analysis")
Standard <- read.csv("SK.csv", header=T, check.names = FALSE)
run1_CT <- read.csv("run1.csv",header=T,check.names = FALSE)
run2_CT <- read.csv("run2.csv",header=T,check.names = FALSE)
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
## run1 ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cal_run1 <- calData(nrmData(data=run1_CT,
                         r=3,
                         E=Ef,
                         Eerror=c(0.02,0.02,0.02,0.02,0.02,0.02), 
                         nSpl=21,
                         nbRef=2,
                         Refposcol=(1:2), 
                         nCTL=5,
                         CF=c(1,1,1,1,1,1), 
                         CalPos=21,
                         trace=FALSE,
                         geo=TRUE,
                         na.rm=TRUE)[[3]])
cal_run1

run1_nrm <- nrmData(data=run1_CT,
                    r=3,
                    E=Ef,
                    Eerror=c(0.02,0.02,0.02,0.02,0.02,0.02), 
                    nSpl=21,
                    nbRef=2,
                    Refposcol=(1:2), 
                    nCTL=5,
                     CF=cal_run1,  
                     CalPos=21, 
                     trace=FALSE,
                     geo=TRUE,
                     na.rm=TRUE)
run1_nrm
run1_nrm_norm <- run1_nrm$`NRQs`
run1_nrm_norm <- run1_nrm_norm[-c(21), ]
run1_nrm_norm
run1_nrm_norm$subspecies <- as.factor(c(rep("carnica", 10),
                                           rep("iberiensis", 10)
))
run1_nrm_norm$role <- as.factor(c(rep("forager bee", 5),
                                            rep("nurse bee", 5),
                                            rep("forager bee", 5),
                                            rep("nurse bee", 5)
))
run1_nrm_norm
run1_nrm_norm_m <- melt(subset(run1_nrm_norm, select=c("subspecies","role","Des1","Des2","El1","El2")), id=c("subspecies","role"), variable.name = "gene", value.name = "expression")
run1_nrm_norm_m

ggplot(data = run1_nrm_norm_m, aes(x=subspecies,
                                   y=expression,
                                   fill=role
                                   )) +
  geom_boxplot(show.legend = TRUE) +
  facet_wrap(~gene,
             scales = "free"
             )+
  theme_classic()+
  labs(y='relative expression')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## run2 ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cal_run2 <- calData(nrmData(data=run2_CT,
                            r=3,
                            E=Ef,
                            Eerror=c(0.02,0.02,0.02,0.02,0.02,0.02), 
                            nSpl=21,
                            nbRef=2,
                            Refposcol=(1:2), 
                            nCTL=5,
                            CF=c(1,1,1,1,1,1), 
                            CalPos=21,
                            trace=FALSE,
                            geo=TRUE,
                            na.rm=TRUE)[[3]])
cal_run2

run2_nrm <- nrmData(data=run2_CT,
                    r=3,
                    E=Ef,
                    Eerror=c(0.02,0.02,0.02,0.02,0.02,0.02), 
                    nSpl=21,
                    nbRef=2,
                    Refposcol=(1:2), 
                    nCTL=5,
                    CF=cal_run2,  
                    CalPos=21, 
                    trace=FALSE,
                    geo=TRUE,
                    na.rm=TRUE)
run2_nrm
run2_nrm_norm <- run2_nrm$`NRQs`
run2_nrm_norm <- run2_nrm_norm[-c(21), ]
run2_nrm_norm
run2_nrm_norm$subspecies <- as.factor(c(rep("carnica", 10),
                                        rep("iberiensis", 10)
))
run2_nrm_norm$role <- as.factor(c(rep("forager bee", 5),
                                  rep("nurse bee", 5),
                                  rep("forager bee", 5),
                                  rep("nurse bee", 5)
))
run2_nrm_norm
run2_nrm_norm_m <- melt(subset(run2_nrm_norm, select=c("subspecies","role","Des1","Des2","El1","El2")), id=c("subspecies","role"), variable.name = "gene", value.name = "expression")
run2_nrm_norm_m

ggplot(data = run2_nrm_norm_m, aes(x=subspecies,
                                   y=expression,
                                   fill=role
)) +
  geom_boxplot(show.legend = TRUE) +
  facet_wrap(~gene,
             scales = "free"
  )+
  theme_classic()+
  labs(y='relative expression')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

