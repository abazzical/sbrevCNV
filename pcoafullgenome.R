
library(colorspace)
library(dbplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library (ape)
library(cowplot)
library(grid)

#
#load necessary packages
AllgenesAllCNVs<- read.table(file="AllgenesAllCNVs.txt", sep = "\t")
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
plotAllGenes<-plot(mountainplotgenes, pch=1, cex = 2 , xlab = "Variance explained = 14.8%", ylab = "Variance explained = 11.0%" , cex.axis = 1.5 , cex.lab = 1.5,xlim=c(-6000, 6000), ylim=c(-4000, 2000)) +
#creates a plot using coast points (this uses gplots so its easier to change cosmetic settings)
points(beachplotgenes, pch = 16 , cex = 2)
#adds mountain points to plot

#to make pcoa w/ sample names:
biplot.pcoa(pcoaall)



# =================================================
#
### make a distance matrix to use in splitstree
### use whole genome to compare with the SNP analysis from Branco et al. (2015) who used all the SNPs

trimDBSplits <- (AllgenesAllCNVs[13:39])
#removes excess data (protein ID's, positions, etc)
transtrimDBSplits <- t(trimDBSplits)
#transposes so samples are y axis
rownames(transtrimDBSplits) <- c("Sb1", "Sb10", "Sb11", "Sb12", "Sb13", "Sb14", "Sb15", "Sb16", "Sb18", "Sb2" ,"Sb20", "Sb21", "Sb22", "Sb23", "Sb24", "Sb25", "Sb26", "Sb27", "Sb28", "Sb29", "Sb3", "Sb34", "Sb4", "Sb5", "Sb6", "Sb8", "Sb9")
#Labels rows by sample
distanceallSplits <- dist(transtrimDBSplits)
#Creates distance matrix from dataset (uses euclidean distance)
forSplits <- as.matrix(distanceallSplits)
write.table(forSplits, "distSplits.txt")
hc = hclust(as.dist(distanceallSplits), method="complete")
phylotree = as.phylo(hc)
plot(phylotree)
write.tree(phy=phylotree, file="tree.nwk") #import this newick tree in splitstree

### calculate pcoa of only regions with at least one absence (a zero)

### filter based on at least one sample has a '0' as estimated number of copies
allSampLoss<-AllgenesAllCNVs %>%
  filter(V39 =="0"|V38=="0"|V37=="0"|V36=="0"|V35=="0"|V34=="0"|V33=="0"|V32=="0"|V31=="0"|V30=="0"|V29=="0"|
           V28=="0"|V27=="0"|V26=="0"|V25=="0"|V24=="0"|V23=="0"|V22=="0"|V21=="0"|V20=="0"|V19=="0"|V18=="0"|
           V17=="0"|V16=="0"|V15=="0"|V14=="0"|V13=="0")
filterallSampLoss <- filter(allSampLoss, V3 == "CDS")
#removes anything that is not a CDS
trimdatabaseLoss <- (filterallSampLoss[13:39])
#removes excess data (protein ID's, positions, etc)
transtrimdatabaseLoss <- t(trimdatabaseLoss)
#transposes so samples are y axis
rownames(transtrimdatabaseLoss) <- c("Sb1", "Sb10", "Sb11", "Sb12", "Sb13", "Sb14", "Sb15", "Sb16", "Sb18", "Sb2" ,"Sb20", "Sb21", "Sb22", "Sb23", "Sb24", "Sb25", "Sb26", "Sb27", "Sb28", "Sb29", "Sb3", "Sb34", "Sb4", "Sb5", "Sb6", "Sb8", "Sb9")
#Labels rows by sample
distanceallLoss <- dist(transtrimdatabaseLoss)#Creates distance matrix from dataset (uses euclidean distance)
pcoaallLoss <- pcoa(distanceallLoss)
#performs PCoA on distance matrix
axesallLoss<- as.matrix(pcoaallLoss[["vectors"]])
axesallLoss <- as.matrix(axesallLoss[,1:2])
#extracts PCA axis 1 and 2 into a data set
mountainplotallLossa <- as.matrix(axesallLoss[1:10,])
mountainplotallLossb <- as.matrix(axesallLoss[23:27,])
mountainplotallLossc <- as.matrix(axesallLoss[21,])
mountainplotallLossc <- t(mountainplotallLossc)
#splits data set into mountains and beaches (extracts mountains)
beachplotallLossa <- as.matrix(axesallLoss[11:20,])
beachplotallLossb <- as.matrix(axesallLoss[22,])
beachplotallLossb <- t(beachplotallLossb)
#extracts beaches
mountainplotgenesLoss <- rbind(mountainplotallLossa, mountainplotallLossb, mountainplotallLossc)
beachplotgenesLoss <- rbind(beachplotallLossa, beachplotallLossb)
#binds datasets into 1 mountain dataset and 1 beach dataset
plotLossGenes<-plot(mountainplotgenesLoss, pch=1, cex = 2 , xlab = "Variance explained = 14.8%", ylab = "Variance explained = 11.0%" , cex.axis = 1.5 , cex.lab = 1.5,xlim=c(-6000, 6000), ylim=c(-2000, 4000))+
#creates a plot using coast points (this uses gplots so its easier to change cosmetic settings)
points(beachplotgenesLoss, pch = 16 , cex = 2)
#adds mountain points to plot

#to make pcoa w/ sample names:
biplot.pcoa(pcoaallLoss)


### calculate pcoa of only regions with at least one gain (more than 2)

### filter based on at least one sample has an estimated number of copies as more than 2
allSampGain<-AllgenesAllCNVs %>%
  filter(V39 >2|V38>2|V37>2|V36>2|V35>2|V34>2|V33>2|V32>2|V31>2|V30>2|V29>2|
           V28>2|V27>2|V26>2|V25>2|V24>2|V23>2|V22>2|V21>2|V20>2|V19>2|V18>2|
           V17>2|V16>2|V15>2|V14>2|V13>2)

filterallSampGain <- filter(allSampGain, V3 == "CDS")
#removes anything that is not a CDS
trimdatabaseGain <- (filterallSampGain[13:39])
#removes excess data (protein ID's, positions, etc)
transtrimdatabaseGain <- t(trimdatabaseGain)
#transposes so samples are y axis
rownames(transtrimdatabaseGain) <- c("Sb1", "Sb10", "Sb11", "Sb12", "Sb13", "Sb14", "Sb15", "Sb16", "Sb18", "Sb2" ,"Sb20", "Sb21", "Sb22", "Sb23", "Sb24", "Sb25", "Sb26", "Sb27", "Sb28", "Sb29", "Sb3", "Sb34", "Sb4", "Sb5", "Sb6", "Sb8", "Sb9")
#Labels rows by sample
distanceallGain <- dist(transtrimdatabaseGain)#Creates distance matrix from dataset (uses euclidean distance)
pcoaallGain <- pcoa(distanceallGain)
#performs PCoA on distance matrix
axesallGain<- as.matrix(pcoaallGain[["vectors"]])
axesallGain <- as.matrix(axesallGain[,1:2])
#extracts PCA axis 1 and 2 into a data set
mountainplotallGaina <- as.matrix(axesallGain[1:10,])
mountainplotallGainb <- as.matrix(axesallGain[23:27,])
mountainplotallGainc <- as.matrix(axesallGain[21,])
mountainplotallGainc <- t(mountainplotallGainc)
#splits data set into mountains and beaches (extracts mountains)
beachplotallGaina <- as.matrix(axesallGain[11:20,])
beachplotallGainb <- as.matrix(axesallGain[22,])
beachplotallGainb <- t(beachplotallGainb)
#extracts beaches
mountainplotgenesGain <- rbind(mountainplotallGaina, mountainplotallGainb, mountainplotallGainc)
beachplotgenesGain <- rbind(beachplotallGaina, beachplotallGainb)
#binds datasets into 1 mountain dataset and 1 beach dataset
plotGainGenes<-plot(mountainplotgenesGain, pch=1, cex = 2 , xlab = "Variance explained = 14.8%", ylab = "Variance explained = 11.0%" , cex.axis = 1.5 , cex.lab = 1.5,xlim=c(-6000, 6000), ylim=c(-4000, 2000))+
#creates a plot using coast points (this uses gplots so its easier to change cosmetic settings)
points(beachplotgenesGain, pch = 16 , cex = 2)
#adds mountain points to plot

#to make pcoa w/ sample names:
biplot.pcoa(pcoaallGain)

