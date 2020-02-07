library(colorspace)
library(dbplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library (ape)
#load necessary packages
filterdatabase1 <- filter(AllgenesAllCNVs, V3 == "CDS")
#removes anything that is not a CDS
trimdatabase1 <- (filterdatabase1[13:39])
#removes excess data (protein ID's, positions, etc)
transposeddatabase <- t(trimdatabase1)
#transposes so samples are y axis
rownames(transposeddatabase) <- c("Sb1", "Sb10", "Sb11", "Sb12", "Sb13", "Sb14", "Sb15", "Sb16", "Sb18", "Sb2" ,"Sb20", "Sb21", "Sb22", "Sb23", "Sb24", "Sb25", "Sb26", "Sb27", "Sb28", "Sb29", "Sb3", "Sb34", "Sb4", "Sb5", "Sb6", "Sb8", "Sb9")
#Labels rows by sample
distanceall <- dist(transposeddatabase)
#Creates distance matrix from dataset (uses euclidean distance)
pcoaall <- pcoa(distanceall)
#performs PCoA on distance matrix
axesall<- as.matrix(pcoaall[["vectors"]])
axesall <- as.matrix(axesall[,1:2])
#extracts PCA axis 1 and 2 into a data set
mountainplotalla <- as.matrix(axesall[1:10,])
mountainplotallb <- as.matrix(axesall[23:27,])
mountainplotallc <- as.matrix(axesall[21,])
mountainplotallc <- t(mountainplotallc)
#splits data set into mountains and beaches (extracts mountains)
beachplotalla <- as.matrix(axesall[11:20,])
beachplotallb <- as.matrix(axesall[22,])
beachplotallb <- t(beachplotallb)
#extracts beaches
mountainplotgenes <- rbind(mountainplotalla, mountainplotallb, mountainplotallc)
beachplotgenes <- rbind(beachplotalla, beachplotallb)
#binds datasets into 1 mountain dataset and 1 beach dataset
plot(mountainplotgenes, pch=1, cex = 2 , xlab = "Variance explained = 14.8%", ylab = "Variance explained = 11.0%" , cex.axis = 1.5 , cex.lab = 1.5)
#creates a plot using coast points (this uses gplots so its easier to change cosmetic settings)
points(beachplotgenes, pch = 16 , cex = 2)
#adds mountain points to plot

#to make pcoa w/ sample names:
biplot.pcoa(pcoaall)
