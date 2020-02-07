library(tidyverse) #dplyr data wrangling
library(stats) #variance function
library(readr) #readcsv
library(cn.mops) #manhattan plot
library(qqman) #manhattan plot
library(graphics) #manhattan plot


# Set working directory to where Ratio files live
  setwd("C:/Users/malth/Documents/CNV/Ratios")
  
# Pathway to folder where Ratio files live 
  folder <- "C:/Users/malth/Documents/CNV/Ratios"
  
#List Ratio files 
  CN_list <-  list.files(path = folder, recursive = TRUE, full.names = FALSE)
 
  
# Read each Ratio file in list 
  CN_read <-  lapply(CN_list, read.table, header=T, sep="\t")

  
  
##TOTAL DATAFRAME 
  
# Create dataframe
  CN_df <- bind_cols(CN_read) %>% #bind ratio files by column
            select(Chromosome, Start, starts_with("Sb"))%>% #select columns: Chromosome, Start, and all Sb Samples 
            as.data.frame() %>% #create dataframe        
            unite("window", Chromosome:Start, remove = T) #create row id
  
  write.csv(CN_df, file = 'C:/Users/malth/Documents/CNV/DataFrame/CN_df.csv') #save dataframe
  
  
# Transpose dataframe
  window <- CN_df$window #remember row id
  CN_total <- as.data.frame(t(CN_df[,-1])) #transpose dataframe, exclude 1st colunmn
  colnames(CN_total) <- window #create column names  
  
write.csv(CN_total, file = 'C:/Users/malth/Documents/CNV/DataFrame/CN_total.csv') #save dataframe

##MONTANE DATAFRAME  
# Create dataframe
  mnt <- c('Sb1', 'Sb2', 'Sb3', 'Sb4', 'Sb5', 'Sb6', 'Sb8', 'Sb9', 'Sb10', 'Sb11', 'Sb12', 'Sb13', 'Sb14', 'Sb15', 'Sb16', 'Sb18')
  CN_mnt <- select(CN_df, window, one_of(mnt))
  
# Transpose dataframe
  window <- CN_mnt$window #remember row id
  CN_mnt <- as.data.frame(t(CN_mnt[,-1])) #transpose dataframe, exclude 1st colunmn
  colnames(CN_mnt) <- window #create column names  
  
write.csv(CN_mnt, file = 'C:/Users/malth/Documents/CNV/DataFrame/CN_mnt.csv') #save dataframe
  
##COASTAL DATAFRAME 
# Create dataframe
  cst <- c('Sb20', 'Sb21', 'Sb22', 'Sb23', 'Sb24', 'Sb25', 'Sb26', 'Sb27', 'Sb28', 'Sb29', 'Sb34')
  CN_cst <- select(CN_df, window, one_of(cst))
  
# Transpose dataframe
  window <- CN_cst$window #remember row id
  CN_cst <- as.data.frame(t(CN_cst[,-1])) #transpose dataframe, exclude 1st colunmn
  colnames(CN_cst) <- window #create column names  
  
write.csv(CN_cst, file = 'C:/Users/malth/Documents/CNV/DataFrame/CN_cst.csv') #save dataframe
  
  
  

##CALCULATE VARIANCE  (Time to run: 2.5hrs / group)
# total variance
  CN_total <- as.tbl(CN_total) #convert dataframe to tibble 
  var_total <- summarise_at(CN_total, vars(1:208889), var) #apply the variance function to each column, output new data frame
  var_total <- t(var_total) #transpose dataframe, output matrix
  colnames(var_total) <- c("Total Variance") #name variable column
  
write.csv(var_total, file = 'C:/Users/malth/Documents/CNV/Variances/var_total.csv') #save dataframe
  
  
# montane variance
  CN_mnt <- as.tbl(CN_mnt) #convert dataframe to tibble 
  var_mnt <- summarise_at(CN_mnt, vars(1:208889), var) #apply the variance function to each column, output new data frame
  var_mnt <- t(var_mnt) #transpose dataframe, output matrix
  colnames(var_mnt) <- c("Total Variance") #name variable column

write.csv(var_mnt, file = 'C:/Users/malth/Documents/CNV/Variances/var_mnt.csv') #save dataframe
  
  
# coastal variance
  CN_cst <- as.tbl(CN_cst) #convert dataframe to tibble 
  var_cst <- summarise_at(CN_cst, vars(1:208889), var) #apply the variance function to each column, output new data frame
  var_cst <- t(var_cst) #transpose dataframe, output matrix
  colnames(var_cst) <- c("Total Variance") #name variable column
  
write.csv(var_cst, file = 'C:/Users/malth/Documents/CNV/Variances/var_cst.csv') #save dataframe
  

  
  
  
##COMBINE VARIANCE DATAFRAMES

# Pathway to folder where variance files live
  variances <- "C:/Users/malth/Documents/CNV/Variances"
  
# List variance files 
  var_list <-  list.files(path = variances, recursive = TRUE, full.names = T)
  
# Read each variance file in the list
  var_read <-  lapply(var_list, read.csv, header=T, sep=",")
  
# Make variances dataframe
  var_df <- bind_cols(var_read) %>% #bind variance files by column
  select(Window, ends_with("Variance"))%>% #select columns to be used in Vst calculation
  as.data.frame()#create dataframe
  
write.csv(var_df, file = 'C:/Users/malth/Documents/CNV/DataFrame/var_df.csv') #save dataframe
  
    
# Calaculate VST
  SbVST <- var_df %>%
          mutate(VST = ((Total.Variance-((Cst.Variance*11+Mnt.Variance*16)/27))/Total.Variance) ) %>% #group.Variance is the name of columns in the variance dataframe
          select(Window, VST) #select columns to be used in Manhattan plot 
  
write.csv(SbVST, file = 'C:/Users/malth/Documents/CNV/DataFrame/SbVST.csv') #save dataframe

  

##In text editor, clean SbVst File 
  # columns : scaffold, position, Vst
  # tabs between columns 
  
  

##MANHATTAN PLOT: Plot all scaffolds, highlight top 1% Vst
#vst

manhattan <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", group = "group", geneNames ="label",
                      col=c("gray10", "gray60"), chrlabs=NULL,
                      suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), 
                      highlight=NULL, highlight1=NULL, highlight2=NULL, highlight3=NULL, highlight4=NULL, highlight5=NULL, 
                      highlight6=NULL, highlight7=NULL, logp=TRUE, annotatePval = NULL, 
                      annotateTop = TRUE, annotate2 = NULL, ...) {
  
  # Not sure why, but package check will warn without this.
  CHR=BP=P=index=NULL
  
  # Check for sensible dataset
  ## Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  ## warn if you don't have a snp column
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # Create a new data.frame with columns called CHR, BP, and P.
  # d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos = NA, index = NA) # with millions of SNPs, create dataframe at once 
  #  rather than dynamically allocated(see line 72-73, and remove line 87 and line 91 )
  
  # If the input data frame has a SNP column, add it to the new data frame you're creating.
  if (!is.null(x[[snp]])) d = data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos = NA, index = NA ,SNP=x[[snp]], stringsAsFactors = FALSE) else 
    d = data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos = NA, index = NA)
  
  
  # Set positions, ticks, and labels for plotting
  ## Sort and keep only values where is numeric.
  #d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
  #  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))       ## unused, all three variables are numeric, line:63-65 
  d <- d[order(d$CHR, d$BP), ]
  #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }
  # d$pos=NA
  
  
  # Fixes the bug where one chromosome is missing by adding a sequential index column.
  # d$index=NA
  # ind = 0
  # for (i in unique(d$CHR)){
  #     ind = ind + 1
  #     d[d$CHR==i,]$index = ind
  # }
  d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP,d$CHR,length))  # replcace the for loop of line 92-96 to improve efficiency
  
  # This section sets up positions and ticks. Ticks should be placed in the
  # middle of a chromosome. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromsome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    ## Uncomment the next two linex to plot single chr results in Mb
    #options(scipen=999)
    #d$pos=d$BP/1e6
    d$pos=d$BP
    #  ticks=floor(length(d$pos))/2+1          ## unused, from code line: 169
    xlabel = paste('Chromosome',unique(d$CHR),'position')
    #  labs = ticks          ## unused, from code line: 169
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced. 
        lastbase = lastbase +max(d[d$index==(i-1),"BP"])   # replace line 128
        d[d$index == i,"BP"] = d[d$index == i,"BP"]-min(d[d$index==i,"BP"]) +1
        d[d$index == i, "pos"] = d[d$index == i,"BP"] + lastbase    # replace line 129
        # lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
        # d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
        
      }
      # Old way: assumes SNPs evenly distributed
      # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
      # New way: doesn't make that assumption
      # ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)  # see line 136, to reduce the burden of for loop 
    }
    ticks <-tapply(d$pos,d$index,quantile,probs=0.5)   # replace line 135
    xlabel = 'Chromosome'
    #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
    labs <- unique(d$CHR)
  }
  
  # Initialize plot
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  
  # The old way to initialize the plot
  # plot(NULL, xaxt='n', bty='n', xaxs='i', yaxs='i', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
  #      xlab=xlabel, ylab=expression(-log[10](italic(p))), las=1, pch=20, ...)
  
  
  # The new way to initialize the plot.=============================================== CHANGE Y-AXIS HERE, DIPSHIT
  ## See http://stackoverflow.com/q/23922130/654296
  ## First, define your default arguments
  def_args <- list(xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                   xlim=c(xmin,xmax), ylim=c(0,0.5),
                   xlab=xlabel, ylab=expression(-log[10](italic(p))))
  ## Next, get a list of ... arguments
  #dotargs <- as.list(match.call())[-1L]
  dotargs <- list(...)
  ## And call the plot function passing NA, your ... arguments, and the default
  ## arguments that were not defined in the ... arguments.
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  
  # If manually specifying chromosome labels, ensure a character vector and number of labels matches number chrs.
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs)==length(labs)) {
        labs <- chrlabs
      } else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    } else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  
  # Add an axis. 
  if (nchr==1) { #If single chromosome, ticks and labels automatic.
    axis(1, ...)
  } else { # if multiple chrs, use the ticks and labels you created above.
    axis(1, at=ticks, labels=labs, ...)
  }
  
  # Create a vector of alternatiting colors
  #col=rep(col, max(d$CHR))  # replaced by line 187
  col = rep_len(col, max(d$index))  ## mean this one?  the results are same
  
  # Add points to the plot
  if (nchr==1) {
    with(d, points(pos, logp, pch=20, col=col[1], ...))
  } else {
    # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
    icol=1
    for (i in unique(d$index)) {
      #with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=20, ...))
      points(d[d$index==i,"pos"], d[d$index==i,"logp"], col=col[icol], pch=20, ...)
      icol=icol+1
    }
  }
  
  # Add suggestive and genomewide lines
  if (suggestiveline) abline(h=suggestiveline, col="blue")
  if (genomewideline) abline(h=genomewideline, col="red")
  ####### CHANGE COLOR HERE ######
  # Highlight snps from a character vector === highlight does the 1% in light blue
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight=d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col="#66B2FF", pch=20, ...)) ##CHANGE COLOR HERE##
  }
  # Highlight ANOTHER group from a character vector ====== YESSSSSSSS highlight2 does 'E-other' category in pale blue
  if (!is.null(highlight1)) {
    if (any(!(highlight1 %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight1=d[which(d$SNP %in% highlight1), ]
    with(d.highlight1, points(pos, logp, col="#3B9AB2", pch=19, cex=1.75, ...)) 
  }
  
  # Highlight ANOTHER group from a character vector ====== YESSSSSSSS highlight2 does 'E-other' category in pale blue
  if (!is.null(highlight2)) {
    if (any(!(highlight2 %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight2=d[which(d$SNP %in% highlight2), ]
    with(d.highlight2, points(pos, logp, col="#3B9AB2", pch=19, cex=1.75, ...)) 
  }
  
  # Highlight ANOTHER group from a character vector === highlight3 does 'A-secreted chelating' 
  if (!is.null(highlight3)) {
    if (any(!(highlight3 %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight3=d[which(d$SNP %in% highlight3), ]
    with(d.highlight3, points(pos, logp, col="#ABDDDE", pch=19, cex=1.75, ...)) 
  }
  # Highlight ANOTHER group from a character vector === highlight4 does 'B-transporters'
  if (!is.null(highlight4)) {
    if (any(!(highlight4 %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight4=d[which(d$SNP %in% highlight4), ]
    with(d.highlight4, points(pos, logp, col="#0B775E", pch=19, cex=1.75, ...)) 
  }
  # Highlight ANOTHER group from a character vector === highlight5 does 'C-chelating inside'
  if (!is.null(highlight5)) {
    if (any(!(highlight5 %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight5=d[which(d$SNP %in% highlight5), ]
    with(d.highlight5, points(pos, logp, col="#35274A", pch=19, cex=1.75, ...)) 
  }
  # Highlight ANOTHER group from a character vector === highlight6 does 'D-Antioxidants'
  if (!is.null(highlight6)) {
    if (any(!(highlight6 %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight6=d[which(d$SNP %in% highlight6), ]
    with(d.highlight6, points(pos, logp, col="#046C9A", pch=19, cex=1.75, ...)) 
  }
  # Highlight ANOTHER group from a character vector === highlight7 does zinc (ZISSOU!) gene SlZnT2 in red, like the beanies
  if (!is.null(highlight7)) {
    if (any(!(highlight7 %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight7=d[which(d$SNP %in% highlight7), ]
    with(d.highlight7, points(pos, logp, col="#F21A00", pch=19, cex=1.75, ...)) 
  }
  
  
  # Annotate subset =============== THIS DONT WORK, I wish it would say: 'you can't use the skeleton arm with that', it would be just as much of an incomprehensible error message, but it would cheer me up.
  if (!is.null(annotate2)) {
    if (any(!(annotate2 %in% d$group))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.annotate2=d[which(d$group %in% annotate2), ]
    with(d.annotate2, textxy(pos, P, offset = 1, labs = d$labels, cex = 0.45)) 
  }
  
  # Highlight top SNPs
  if (!is.null(annotatePval)) {
    # extract top SNPs at given p-val
    if (logp) {
      topHits = subset(d, P <= annotatePval)
    } else
      topHits = subset(d, P >= annotatePval)
    par(xpd = TRUE)
    # annotate these SNPs
    if (annotateTop == FALSE) {
      if (logp) {
        with(subset(d, P <= annotatePval), 
             textxy(pos, -log10(P), offset = 0.625, labs = topHits$SNP, cex = 0.45), ...)
      } else
        with(subset(d, P >= annotatePval), 
             textxy(pos, P, offset = 0.625, labs = topHits$SNP, cex = 0.45), ...)
    }
    else {
      # could try alternative, annotate top SNP of each sig chr
      topHits <- topHits[order(topHits$P),]
      topSNPs <- NULL
      
      for (i in unique(topHits$CHR)) {
        
        chrSNPs <- topHits[topHits$CHR == i,]
        topSNPs <- rbind(topSNPs, chrSNPs[1,])
        
      }
      if (logp ){
        textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, labs = topSNPs$SNP, cex = 0.5, ...)
      } else
        textxy(topSNPs$pos, topSNPs$P, offset = 0.625, labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }  
  par(xpd = FALSE)
}


# load VST Dataframe  
  vst <- read.csv("C:/Users/malth/Documents/CNV/DataFrames/SbVST.csv", sep="\t")
  
# set up dataframe for manhattan plot function
  vst <- filter(vst, Vst != "NA") #remove rows where Vst values is NA
  vst$scaffold <- as.numeric(vst$scaffold) #change dataclass of scaffold column to numeric
  vst$Vst <- as.numeric(vst$Vst) #change dataclass of Vst column to numeric 
  vst_final<- mutate(vst, group = case_when(Vst > 0.1746032 ~ "grp1", Vst <0.1746032 ~ "other")) #identify the top 1% of Vst values, but assigning each window to a group
  pct1=c("grp1") #define the 1% group 
  
# Plot total + highlight top 1%
  manhattan(vst_final, chr = "scaffold", bp = "position", p = "Vst", snp = "group", logp=FALSE, 
            ylab="Vst", xlab = "scaffold", col = c("gray10","gray60"), ylim =  c(0, 1), cex.axis = 0.9,
            highlight = pct1)
  
# Plot scaffolds 1-100 + highlight top 1%
  vst100 <- filter(vst, scaffold <= 100, Vst != "NA")
  vst100_final <- mutate(vst100, group = case_when(Vst > 0.10238095 ~ "grp1", Vst <0.10238095 ~ "other")) #identify the 1% values, by assigning each scaffold window to a group
  pct1=c("grp1") #define the 1% group
  
  manhattan(vst100_final, chr = "scaffold", bp = "position", p = "Vst", snp = "group", logp=FALSE, ylab="Vst", xlab = "scaffold",
            col = c("gray10","gray60"), ylim =  c(0, 0.4), cex.axis = 0.9,
            highlight = pct1)
  
## scaffolds 1-100 proportion of the genome
  
# Length of whole genome  
  vst_n <- vst %>%
    group_by(scaffold)%>%
    count()
  vst_n_bp <- sum(vst_n$n)*250
  print(vst_n_bp)
  
# Length of scaffolds 1-100  
  vst100_n <- filter(vst, scaffold <= 100)%>%
    group_by(scaffold)%>%
    count()
  vst100_n_bp <- sum(vst100_n$n)*250
  print(vst100_n_bp)
  
# Proportion of genome - scaffolds 1-100
  vst100_bp_ratio<- (vst100_n_bp/vst_n_bp)
  print(vst100_bp_ratio)
  
  



##EXTRA MANHATTAN PLOT CODE 

## Bar Plot of Average Vst/scaffold
vst_avg <- vst %>%
  group_by(scaffold)%>% 
  summarise(mean(Vst))

barplot(vst_avg$`mean(Vst)`, ylim = c(0,0.25),xpd = F, xlab = "scaffold", ylab = "Average Vst") #plot 
abline(v=101, col = "red") #add line to show where scaffolds 1-100 are located 


## dataframe of top 1% Vst values of all scaffolds
vst1<-filter(vst, Vst > 0.1746032)

write.csv(vst1, file = 'C:/Users/malth/Documents/CNV/DataFrame/vst1.csv') #save dataframe


## dataframe of top 1% Vst across all scaffolds, to be used for heat map
vst4<-filter(vst, Vst > 0.1746032) %>%
  mutate(end = position + 250)%>%
  rename(position = "start")%>%
  select(scaffold, start, end, Vst)

write.csv(vst4, file = 'C:/Users/malth/Documents/CNV/vst4.csv') #save dataframe


## dataframe of scaffolds 1-100 
vst100 <- filter(vst, scaffold <= 100 )

write.csv(vst, file = 'C:/Users/malth/Documents/CNV/DataFrame/vst100.csv') #save dataframe


## dataframe of the top 1% Vst values for scaffolds 1-100
vst100_1 <- filter(vst100, Vst > 0.10238095)

write.csv(vst100_1, file = 'C:/Users/malth/Documents/CNV/DataFrame/vst100_1.csv') #save dataframe


## Plot scaffolds with length >=3000bp
# determine which scaffolds have more than 3000bp
vst_gp3 <- filter(vst, position > 30000, Vst != "NA")%>% 
  group_by(scaffold)%>%
  count() 

# pull out scaffolds that have more than 3000bp
vst3000 <- vst %>% 
  semi_join(vst_gp3)%>%
  filter(Vst != "NA")

vst3000$scaffold <- as.numeric(vst3000$scaffold)
vst3000$Vst <- as.numeric(vst3000$Vst)

manhattan(vst3000, chr = "scaffold", bp = "position", p = "Vst", logp=FALSE, ylab="Vst", xlab = "scaffold",
          col = c("gray10","gray60"), ylim =  c(0, 0.4), cex.axis = 0.9)




