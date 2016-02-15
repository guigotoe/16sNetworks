####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: February 2016
# Created: February 2016
#
# This is written as part of 16S - Aging analysis, but could be
# splitted to serve different purposes.
####################################################

#* requirements *#
library(ggplot2)
library(reshape)
library(reshape2)
library(scales)
library(vegan)
library(ade4)
library(RColorBrewer) 

script_loc = getwd()
output_loc = "/Users/guillermotorres/Documents/Proyectos/Doctorado/Metagenome_16S/16sNetworks/results/plotsMothur/"
files.path = '/Users/guillermotorres/Documents/Proyectos/Doctorado/Metagenome_16S/16sNetworks/results/fromMothur/'#
setwd(output_loc)

counts <- read.table(paste(files.path,'16S.an.shared',sep=''),header=T,sep="\t")
taxonomy <- read.table(paste(files.path,'16S.an.0.03.cons.taxonomy',sep=''),header=T,sep="\t")
summary <- read.table(paste(files.path,'16S.an.groups.summary',sep=''),header=T,sep="\t")

#### 
threshold <- quantile(summary$ace,probs=0.85)  ## Estimated population size based on rare OTUs calculated by ace
n=50
otu_size_dist <- data.frame("Size"=character(n),"OTUs"= numeric(n),"more_OTUs"= numeric(n),"less_OTUs"= numeric(n),stringsAsFactors = FALSE)
for (i in seq(1:n)){
  otu_size_dist[i,] <- c(i,table(taxonomy[2]==i)["TRUE"],table(taxonomy[2]>=i)["TRUE"],table(taxonomy[2]<=i)["TRUE"])
}
otu_size_dist[n,]$Size <- paste(n,'+/-',sep=" ")
min_read_len <- tail(which(otu_size_dist$more_OTUs >= threshold),1)
plot(otu_size_dist[,3],main="Distribution of OTU sizes",ylab="Num. OTUs with more than 'x' Num. sequences",xlab="Num. sequences")
abline(v=otu_size_dist$Size[min_read_len],col="red",lty=2)
trim_otus <- otu_size_dist[min_read_len:NROW(otu_size_dist),3]
plot(trim_otus,,main="Distribution of OTU sizes - trimmed data",ylab="Num. OTUs with more than 'x' Num. sequences",xlab="Num. sequences")
sum(taxonomy$Size[taxonomy$Size>=as.integer(otu_size_dist$Size[min_read_len])])/sum(taxonomy$Size) # retained info 0.9667
require(fitdistrplus)
otu_dist <- fitdist(otu_size_dist[,3],"lnorm")
otu_dist
plotdist(otu_size_dist[,3],"lnorm",para=list(meanlog=otu_dist$estimate[1],sdlog=otu_dist$estimate[2]))
otu_trim_dist <- fitdist(trim_otus,"lnorm")
otu_trim_dist
plotdist(trim_otus,"lnorm",para=list(meanlog=otu_trim_dist$estimate[1],sdlog=otu_trim_dist$estimate[2]))
####

####
read_threshold = as.integer(otu_size_dist$Size[min_read_len])
newtax <- taxonomy[taxonomy$Size>=read_threshold,]
newcounts <- counts[,-c(1,2,3)][,newtax$OTU]
rownames(newcounts) <- counts[,2]


