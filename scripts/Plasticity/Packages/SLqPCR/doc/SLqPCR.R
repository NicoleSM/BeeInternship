### R code from vignette source 'SLqPCR.Rnw'

###################################################
### code chunk number 1: SLqPCR
###################################################
library(SLqPCR)
data(vandesompele)
str(vandesompele)


###################################################
### code chunk number 2: fig2
###################################################
tissue <- as.factor(c(rep("BM", 9), rep("POOL", 9), rep("FIB", 20), rep("LEU", 13), rep("NB", 34)))
res.BM <- selectHKgenes(vandesompele[tissue == "BM",], method = "Vandesompele", geneSymbol = names(vandesompele), minNrHK = 2, trace = TRUE, na.rm = TRUE)
res.POOL <- selectHKgenes(vandesompele[tissue == "POOL",], method = "Vandesompele", geneSymbol = names(vandesompele), minNrHK = 2, trace = FALSE, na.rm = TRUE)
res.FIB <- selectHKgenes(vandesompele[tissue == "FIB",], method = "Vandesompele", geneSymbol = names(vandesompele), minNrHK = 2, trace = FALSE, na.rm = TRUE)
res.LEU <- selectHKgenes(vandesompele[tissue == "LEU",], method = "Vandesompele", geneSymbol = names(vandesompele), minNrHK = 2, trace = FALSE, na.rm = TRUE)
res.NB <- selectHKgenes(vandesompele[tissue == "NB",], method = "Vandesompele", geneSymbol = names(vandesompele), minNrHK = 2, trace = FALSE, na.rm = TRUE)


###################################################
### code chunk number 3: table3
###################################################
ranks <- data.frame(c(1, 1:9), res.BM$ranking, res.POOL$ranking, res.FIB$ranking, res.LEU$ranking, res.NB$ranking)
names(ranks) <- c("rank", "BM", "POOL", "FIB", "LEU", "NB")
ranks


###################################################
### code chunk number 4: relQuant
###################################################
exa1 <- apply(vandesompele[tissue == "BM",], 2, relQuantPCR, E = 2)


###################################################
### code chunk number 5: fig2
###################################################
library(RColorBrewer)
mypalette <- brewer.pal(5, "Set1")
matplot(cbind(res.BM$meanM, res.POOL$meanM, res.FIB$meanM, res.LEU$meanM, res.NB$meanM), type = "b", ylab = "Average expression stability M", xlab = "Number of remaining control genes", axes = FALSE, pch = 19, col = mypalette, ylim = c(0.2, 1.22), lty = 1, lwd = 2, main = "Gene stability measure")
axis(1, at = 1:9, labels = as.character(10:2))
axis(2, at = seq(0.2, 1.2, by = 0.2), labels = as.character(seq(0.2, 1.2, by = 0.2)))
box()
abline(h = seq(0.2, 1.2, by = 0.2), lty = 2, lwd = 1, col = "grey")
legend("topright", legend = c("BM", "POOL", "FIB", "LEU", "NB"), fill = mypalette)


###################################################
### code chunk number 6: fig3a
###################################################
mypalette <- brewer.pal(8, "YlGnBu")
barplot(cbind(res.BM$variation, res.POOL$variation, res.FIB$variation, res.LEU$variation, res.NB$variation), beside = TRUE, col = mypalette, space = c(0, 2), names.arg = c("BM", "POOL", "FIB", "LEU", "NB"))
legend("topright", legend = c("V9/10", "V8/9", "V7/8", "V6/7", "V5/6", "V4/5", "V3/4", "V2/3"), fill = mypalette, ncol = 2)
abline(h = seq(0.05, 0.25, by = 0.05), lty = 2, col = "grey")
abline(h = 0.15, lty = 1, col = "black")


###################################################
### code chunk number 7: norm
###################################################
data(SLqPCRdata)
SLqPCRdata
(relData <- apply(SLqPCRdata, 2, relQuantPCR, E = 2))
geneStabM(relData[,c(3,4)])
(exprData <- normPCR(SLqPCRdata, c(3,4)))


