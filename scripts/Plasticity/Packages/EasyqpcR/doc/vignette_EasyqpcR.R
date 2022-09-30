### R code from vignette source 'vignette_EasyqpcR.Rnw'

###################################################
### code chunk number 1: first
###################################################

library(EasyqpcR)

data(Efficiency_calculation)

slope(data=Efficiency_calculation, q=c(1000, 100 ,10, 1, 0.1),
    r=3, na.rm=TRUE)



###################################################
### code chunk number 2: step1
###################################################

efficiency <- slope(data=Efficiency_calculation, q=c(1000, 100 ,10, 1, 0.1),
    r=3, na.rm=TRUE)



###################################################
### code chunk number 3: step2
###################################################
data(qPCR_run1,qPCR_run2,qPCR_run3)

str(c(qPCR_run1,qPCR_run2,qPCR_run3))


###################################################
### code chunk number 4: step3
###################################################

## Isolate the calibrator NRQ values of the first biological replicate

aa <- nrmData(data=qPCR_run1 , r=3, E=c(2, 2, 2, 2),
       Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5,
       nbRef=2, Refposcol=1:2, nCTL=2,
       CF=c(1, 1, 1, 1), CalPos=5, trace=FALSE, geo=TRUE, na.rm=TRUE)[[3]]

## Isolate the calibrator NRQ values of the first biological replicate

bb <- nrmData(data=qPCR_run2 , r=3, E=c(2, 2, 2, 2),
	   Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5,
	   nbRef=2, Refposcol=1:2, nCTL=2,
	   CF=c(1, 1, 1, 1), CalPos=5, trace=FALSE, geo=TRUE, na.rm=TRUE)[[3]]

## Isolate the calibrator NRQ values of the first biological replicate

cc <- nrmData(data=qPCR_run3 , r=3, E=c(2, 2, 2, 2),
       Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5,
       nbRef=2, Refposcol=1:2, nCTL=2,
       CF=c(1, 1, 1, 1), CalPos=5, trace=FALSE, geo=TRUE, na.rm=TRUE)[[3]]



###################################################
### code chunk number 5: step4
###################################################

## Calibration factor calculation

e <- calData(aa)

f <- calData(bb)

g <- calData(cc)



###################################################
### code chunk number 6: step5 (eval = FALSE)
###################################################
## 
## nrmData(data=qPCR_run1 , r=3, E=c(2, 2, 2, 2),
##        Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5,
##        nbRef=2, Refposcol=1:2, nCTL=2,
##        CF=e, CalPos=5, trace=FALSE, geo=TRUE, na.rm=TRUE)
## 
## nrmData(data=qPCR_run2 , r=3, E=c(2, 2, 2, 2),
##        Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5,
##        nbRef=2, Refposcol=1:2, nCTL=2,
##        CF=f, CalPos=5, trace=FALSE, geo=TRUE, na.rm=TRUE)
## 
## nrmData(data=qPCR_run3 , r=3, E=c(2, 2, 2, 2),
##        Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5,
##        nbRef=2, Refposcol=1:2, nCTL=2,
##        CF=g, CalPos=5, trace=FALSE, geo=TRUE, na.rm=TRUE)
## 


###################################################
### code chunk number 7: step6
###################################################

## Isolate the NRQs scaled to control of the first biological replicate

a1 <- nrmData(data=qPCR_run1 , r=3, E=c(2, 2, 2, 2),
       Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5,
       nbRef=2, Refposcol=1:2, nCTL=2,
       CF=e, CalPos=5, trace=FALSE, geo=TRUE, na.rm=TRUE)[1]

## Isolate the NRQs scaled to control of the second biological replicate

b1 <- nrmData(data=qPCR_run2 , r=3, E=c(2, 2, 2, 2),
       Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5,
       nbRef=2, Refposcol=1:2, nCTL=2,
       CF=f, CalPos=5, trace=FALSE, geo=TRUE, na.rm=TRUE)[1]

## Isolate the NRQs scaled to control of the third biological replicate

c1 <- nrmData(data=qPCR_run3 , r=3, E=c(2, 2, 2, 2),
       Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5,
       nbRef=2, Refposcol=1:2, nCTL=2,
       CF=g, CalPos=5, trace=FALSE, geo=TRUE, na.rm=TRUE)[1]

## Data frame transformation

a2 <- as.data.frame(a1)
b2 <- as.data.frame(b1)
c2 <- as.data.frame(c1)

## Aggregation of the three biological replicates

d2 <- rbind(a2, b2, c2)



###################################################
### code chunk number 8: step7
###################################################

totData(data=d2, r=3, geo=TRUE, logarithm=TRUE, base=2,
       transformation=TRUE, nSpl=5, linear=TRUE,
       na.rm=TRUE)



###################################################
### code chunk number 9: step8
###################################################

file <- system.file("extdata", "qPCR_run1.csv", package="EasyqpcR")

qPCR_run1 <- read.table(file, header=TRUE, sep="", dec=".")

qPCR_run1


###################################################
### code chunk number 10: step9
###################################################

badCt(data=qPCR_run1, r=3, threshold=0.5, na.rm=TRUE)



###################################################
### code chunk number 11: step10
###################################################

badCt(data=qPCR_run1, r=3, threshold=0.2, na.rm=TRUE)



###################################################
### code chunk number 12: step11
###################################################

filebis <- system.file("extdata", "Gene_maximisation.csv", package="EasyqpcR")

Gene_maximisation <- read.table(filebis, header=TRUE, sep=";", dec=",")



###################################################
### code chunk number 13: step12
###################################################

badCt(data=Gene_maximisation, r=3, threshold=0.5, na.rm=FALSE)[1]



###################################################
### code chunk number 14: step13 (eval = FALSE)
###################################################
## 
## fileter <- system.file("extdata", "Gene_maximisation_cor.csv",
##     package="EasyqpcR")
## 
## Gene_maximisation_cor <- read.table(fileter, header=TRUE, sep=";", dec=",")
## 
## Gene_maximisation_cor1 <- Gene_maximisation_cor[-c(106:108, 118:120, 130:132,
##  142:144, 154:156, 166:168, 178:180, 190:192),]
## 
## rownames(Gene_maximisation_cor1) <- c(1:168)
## 


###################################################
### code chunk number 15: step14 (eval = FALSE)
###################################################
## 
## calr1 <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2,
##     nCTL=16, CF=c(1, 1, 1), CalPos=c(33:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[3]][1:3,]
## 
## calr2 <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=c(1, 1, 1), CalPos=c(33:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[3]][4:6,]
## 
## calr3 <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=c(1, 1, 1), CalPos=c(33:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[3]][7:9,]
## 
## calr4 <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2,
##     nCTL=16, CF=c(1, 1, 1), CalPos=c(33:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[3]][10:12,]
## 
## calr5 <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=c(1, 1, 1), CalPos=c(33:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[3]][13:15,]
## 
## calr6 <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=c(1, 1, 1), CalPos=c(33:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[3]][16:18,]
## 
## calr7 <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2,
##     nCTL=16, CF=c(1, 1, 1), CalPos=c(33:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[3]][19:21,]
## 
## calr8 <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=c(1, 1, 1), CalPos=c(33:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[3]][22:24,]
## 
## 
## e <- calData(calr1)
## 
## f <- calData(calr2)
## 
## g <- calData(calr3)
## 
## h <- calData(calr4)
## 
## i <- calData(calr5)
## 
## j <- calData(calr6)
## 
## k <- calData(calr7)
## 
## l <- calData(calr8)
## 
## 


###################################################
### code chunk number 16: step15 (eval = FALSE)
###################################################
## 
## m <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=e, CalPos=c(33:35), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[2]][c(1:4,33:35),]
## 
## 
## n <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=f, CalPos=c(36:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[2]][c(5:8,36:38),]
## 
## o <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=g, CalPos=c(36:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[2]][c(9:12,39:41),]
## 
## p <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=h, CalPos=c(33:35), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[2]][c(13:16,42:44),]
## 
## 
## q <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=i, CalPos=c(36:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[2]][c(17:20,45:47),]
## 
## r <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=j, CalPos=c(36:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[2]][c(21:24,48:50),]
## 
## s <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=k, CalPos=c(33:35), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[2]][c(25:28,51:53),]
## 
## 
## t <- nrmData(data = Gene_maximisation_cor1, r=3, E=c(2, 2, 2),
##     Eerror=c(0.02, 0.02, 0.02), nSpl=56, nbRef=2, Refposcol=1:2, nCTL=16,
##     CF=l, CalPos=c(36:56), trace = FALSE, geo = TRUE,
##     na.rm = TRUE)[[2]][c(29:32,54:56),]
## 
## ## Aggregation of all the CNRQs
## 
## u <- rbind(m, n, o, p, q, r, s, t)
## 


###################################################
### code chunk number 17: step16 (eval = FALSE)
###################################################
## 
## ctlgroup <- u[c(1:4,8:11,15:18,22:25),]
## 
## ctlgeom <- colProds(ctlgroup)^(1/dim(ctlgroup)[1])
## ctlgeom1 <- (as.data.frame(ctlgeom)[rep(1:(ncol(u)), each = nrow(u)), ])
## ctlgeom2 <- as.data.frame(matrix(ctlgeom1, ncol = ncol(u), byrow = FALSE))
## 
## CNRQs_scaled_to_group <- u/ctlgeom2
## 
## 


