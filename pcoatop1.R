library(colorspace)
library(dbplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library (ape)
#load necessary packages
dataframe3 <- top1genesforpcoa
dataframe3 <- t(dataframe3)
dataframe3 <- as.data.frame(dataframe3[,2:24018])
#Creates dataframe from copy number variation of top 1% of genes by VST
rownames(dataframe3) <- c("Sb1", "Sb10", "Sb11", "Sb12", "Sb13", "Sb14", "Sb15", "Sb16", "Sb18", "Sb2" ,"Sb20", "Sb21", "Sb22", "Sb23", "Sb24", "Sb25", "Sb26", "Sb27", "Sb28", "Sb29", "Sb3", "Sb34", "Sb4", "Sb5", "Sb6", "Sb8", "Sb9")
#Labels rows by sample
distance3 <- dist(dataframe3)
#Creates distance matrix from dataset (uses euclidean distance)
pcoatop1 <- pcoa(distance3)
#performs PCoA on distance matrix
axestop1<- as.matrix(pcoatop1[["vectors"]])
axestop1 <- as.matrix(axestop1[,1:2])
#extracts PCA axis 1 and 2 into a data set
mountainplotaa <- as.matrix(axestop1[1:10,])
mountainplotbb <- as.matrix(axestop1[23:27,])
mountainplotcc <- as.matrix(axestop1[21,])
mountainplotcc <- t(mountainplotcc)
#splits data set into mountains and beaches (extracts mountains)
beachplotaa <- as.matrix(axestop1[11:20,])
beachplotbb <- as.matrix(axestop1[22,])
beachplotbb <- t(beachplotbb)
#extracts beaches
mountainplottop1 <- rbind(mountainplotaa, mountainplotbb, mountainplotcc)
beachplottop1 <- rbind(beachplotaa, beachplotbb)
#binds datasets into 1 mountain dataset and 1 beach dataset
plot(mountainplottop1, pch=1 , cex = 2 , xlab = "Variance explained = 14.8%", ylab = "Variance explained = 11.0%", cex.axis = 1.5 , cex.lab = 1.5 , cex.point = 2)
#creates a plot using coast points
points(beachplottop1, pch = 16, cex = 2 )
#adds mountain points to plot

#to make pcoa w/ sample names:
biplot.pcoa(pcoatop1)

