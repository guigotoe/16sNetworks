####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: 30 October 2015
# Created: 29 October 2015
#
# This is written as part of 16S - Aging analysis, but could be
# splitted to serve different purposes.
####################################################

#* requirements *#
library(ggplot2)
library(ggrepel)
library(reshape2)
library(scales)
library(vegan)
library(ade4)
library(RColorBrewer) 

source("/home/torres/Documents/Projects/Metagenome/bin/rscripts/16sFunctions.R")

#* Globals *#
script_loc = getwd()
setwd("/home/torres/Documents/Projects/Metagenome/results/plotsMothur/03.2016/")
files.path = '/home/torres/ikmb_storage/Metagenome/16s/03.2016/'#
send.to <- "/home/torres/Documents/Projects/Metagenome/MothurResults/03.2016/"
plots <- "/home/torres/Documents/Projects/Metagenome/results/plotsMothur/03.2016/"




obs_files <- list.files(files.path,pattern="\\.sobs$")
rarefaction <- read.table(paste(files.path,'16S.an.groups.rarefaction',sep=''),header=T,sep="\t")
rarefactionN <- read.table(paste(files.path,'16S.an.0.03.subsample.groups.rarefaction',sep=''),header=T,sep="\t")
taxonomy <- read.table(paste(files.path,'16S.an.0.03.cons.taxonomy',sep=''),header=T,sep="\t")

calcN <- read.table(paste(files.path,'16S.an.0.03.subsample.groups.summary',sep=''),header=T,sep="\t")
otu_ab <- read.table(paste(files.path,'16S.an.cons.taxonomy',sep=''),header=T,sep="\t")
counts <- read.table(paste(files.path,'16S.an.0.03.subsample.shared',sep=''),header=T,sep="\t") 
wdeg <- read.table('wdegree.txt',header=T,sep="\t")
ass <- read.table('assortcoef.txt',header=T,sep="\t")
acc <- read.table('avclustcoef.txt',header=T,sep="\t")
asp <- read.table('avshortest.txt',header=T,sep="\t")
bwc <- read.table('betwcent.txt',header=T,sep="\t")
clos <- read.table('closeness.txt',header=T,sep="\t")
eigc <- read.table('eigvectc.txt',header=T,sep="\t")


#################################
#* diversity                   *#
#################################


#############################
## Colector && Rarefaction ##

#sobs.df <- OTUs(obs_files) # Colector curve using Observed OTUs

## Colector && Rarefaction ##
rarefaction <- read.table(paste(files.path,'16s.an.groups.rarefaction',sep=''),header=T,sep="\t")
rarefaction$X <- NULL
rare.df <- data.frame("numsampled"=rarefaction[1:which(rarefaction[1]==15000),1])  # Rarefaction curve calculated using Observed OTUs
for (i in seq(from=1,to=length(rarefaction[,-1]),by=3)){rare.df <- cbind(rare.df,rarefaction[1:which(rarefaction[1]==15000),-1][i])}
rare.graph(rare.df,by.g=T,lg=T,savef=plots)

##Rarefaction subsample - normalized
rarefactionN <- read.table(paste(files.path,'16S.an.0.03.subsample.groups.rarefaction',sep=''),header=T,sep="\t")
rarefactionN$X <- NULL
raren.df <- rarefactionN[1]  # Rarefaction curve calculated using Observed OTUs
for (i in seq(from=1,to=length(rarefactionN[,-1]),by=3)){raren.df <- cbind(raren.df,rarefactionN[,-1][i])}
rare.graph(raren.df,by.g=T,lg=F)

######
# End #
#######

##Alpha diversity normalized
# invsimpson and Shannon entropy

idxs <- read.table(paste(files.path,'16s.an.groups.summary',sep=''),header=T,sep="\t") # importing summary file with all indexes
design <- read.table(paste(files.path,'design.txt',sep=''),header=T,sep="\t",row.names=1)
idxs$X <- NULL
divs <- add.GenderAge(idxs,design,colname="group")
div.graphs(divs,savef=plots)


############################################
#* Data filtering and Taxonomic Analysis  *#
############################################

source("/home/torres/Documents/Projects/Metagenome/bin/rscripts/16sFunctions.R")
taxonomy <- read.table(paste(files.path,'16s.an.cons.taxonomy',sep=''),header=T,sep="\t")
counts <- read.table(paste(files.path,'16s.an.shared',sep=''),header=T,sep="\t")
counts$X <- NULL
idxs <- read.table(paste(files.path,'16s.an.groups.summary',sep=''),header=T,sep="\t") # importing summary file with all indexes
design <- read.table(paste(files.path,'design.txt',sep=''),header=T,sep="\t",row.names=1)
idxs$X <- NULL

raw_counts <- add.GenderAge(counts[,2:3],design,colname="Group") # retrieve ID,Gender,Age,Age.group info
raw_counts$numOtus <- NULL
raw_counts <- merge(raw_counts,counts,by="Group",all.x=T)

min_bin_len <- otu.filtering(taxonomy,qval=0.85,nval=100,plots=F) ## for more info check the function... update for plots

##################################
## *  New files from filter  * ##
################################
newtax <- taxonomy[taxonomy$Size>=min_bin_len,]
newraw_counts <- raw_counts[,-c(1,7,8)][,newtax$OTU]
rownames(newraw_counts) <- raw_counts[,1]

#################
# * Tax. plots #
###############

rdata <- newraw_counts[,-c(1:5)]
read_retained <- sum(as.numeric(rdata))/sum(as.numeric(unlist(raw_counts[,-c(1:8)]))) # 0.9717457

# get taxa category
rdata["phylla",] <- unlist(lapply(colnames(rdata),function(x) taxa(x,2))) 
rdata["genus",] <- unlist(lapply(colnames(rdata),function(x) taxa(x,6)))
rdata["family",] <- unlist(lapply(colnames(rdata),function(x) taxa(x,5)))

write.table(rdata,paste(send.to,"rdata",sep=""),sep="\t",col.names=T,row.names=T,quote="")
rdata[is.na(rdata)] <- "unclassified"
write.table(rdata,paste(send.to,"rdata_un",sep=""),sep="\t",col.names=T,row.names=T,quote="")

phylla.raw <- rdata
genus.raw <- rdata
family.raw <- rdata

## Removing unclasified taxlevels
for (i in seq(1:length(rdata["phylla",]))){
  if(rdata["phylla",i]=="unclassified") phylla.raw[,colnames(rdata[i])] <- NULL
  if(rdata["family",i]=="unclassified") family.raw[,colnames(rdata[i])] <- NULL
  if(rdata["genus",i]=="unclassified") genus.raw[,colnames(rdata[i])] <- NULL
}
#sum(data)/sum(rdata)
source("/home/torres/Documents/Projects/Metagenome/bin/rscripts/16sFunctions.R")
#################
##** Phylla **##

## without unclassified
phylla.dfs <-  getTaxGenderDF(df.raw=phylla.raw,taxlevel="phylla",design=design)
tax_graph(phylla.dfs$all,taxlevel="Phylla",savef=plots,ids=F)
## with unclassified
phylla.dfs.plus <-  getTaxGenderDF(df.raw=rdata,taxlevel="phylla",design=design)
tax_graph(phylla.dfs.plus$all,taxlevel="Phylla",savef=plots,ids=F)

#############
# * Famlily

## without unclassified
family.dfs <-  getTaxGenderDF(df.raw=family.raw,taxlevel="family",design=design)
tax_graph(family.dfs$all,taxlevel="Families",savef=plots,ids=F)
## with unclassified
family.dfs.plus <-  getTaxGenderDF(df.raw=rdata,taxlevel="family",design=design)
tax_graph(family.dfs.plus$all,taxlevel="Families",savef=plots,ids=F)

###############
##** GENUS **##

## without unclassified
genus.dfs <-  getTaxGenderDF(df.raw=genus.raw,taxlevel="genus",design=design)
tax_graph(genus.dfs$all,taxlevel="Genus",savef=plots,ids=F)
## with unclassified
genus.dfs.plus <-  getTaxGenderDF(df.raw=rdata,taxlevel="genus",design=design)
tax_graph(genus.dfs.plus$all,taxlevel="Genus",savef=plots,ids=F)

####
dim(phylla.raw)
dim(family.raw)
dim(genus.raw)

#####################################################################################
# * Compositional Data treatment; Outliers, strong PCA and Discriminant analysis * ##
#####################################################################################
df.list <- phylla.dfs
taxlevel <- "Phylla"
cpanalysis(df.list,taxlevel,savef=plots){

#################################
#* Networks metrics graphs     *#
#################################

head(wdeg[,1:5])
qplot(x=wdeg$'X',y=wdeg$'Actinomycetaceae_Actinomyces_Otu001125',data=wdeg,geom="line")
rownames(wdeg) = wdeg$X
wdm <- apply(wdeg[2:ncol(wdeg)],1,mean)
head(wdm)
wdeg$mean <- wdm
ggplot(data=wdeg,aes(X,mean))+geom_point()

ggplot(data=acc,aes(X,ave_cluster_coef))+geom_point()
ggplot(data=asp,aes(X,ave_cluster_coef))+geom_point()+ylab("ave_shortest_path")

head(bwc[,1:5])
bwcm <- apply(bwc[2:ncol(bwc)],1,mean)
bwc$mean <- bwcm
ggplot(data=bwc,aes(X,mean))+geom_point()#+ylab("mean_closness")

head(clos[,1:5])
closm <- apply(clos[2:ncol(clos)],1,mean)
clos$mean <- closm
ggplot(data=clos,aes(X,mean))+geom_point()+ylab("mean_closness")

head(eigc[,1:5])
eigcm <- apply(eigc[2:ncol(eigc)],1,mean)
eigc$mean <- eigcm
ggplot(data=eigc,aes(X,mean))+geom_point()+ylab("mean_eigenvector_centrality")









#################
##### Others ####

#################################
#* Reads in OTU - Distribution *#
#################################

source("/home/torres/Documents/Projects/Metagenome/bin/rscripts/16sFunctions.R")
taxonomy <- read.table(paste(files.path,'16S.an.cons.taxonomy',sep=''),header=T,sep="\t")
counts <- read.table(paste(files.path,'16s.an.shared',sep=''),header=T,sep="\t")
counts$X <- NULL
head(counts[,1:5])

df <- add.GenderAge(counts[,2:3],design,colname="Group") # retrieve ID,Gender,Age,Age.group info
df$numOtus <- NULL
df <- merge(df,counts,by="Group",all.x=T)

tdf <- t(df[,8:ncol(df)])
colnames(tdf) <- df$
  head(tcount[,1:5])
tcount <- as.data.frame(tcount)
head(tcount[,1:5])

## calculating OTUs per individuals ##
opi <- data.frame(indiv=integer(),OTUs=integer())
ones <- apply(tcount,2,function(x) x!=0)
otuN <- apply(ones,1,function(x) table(x)["TRUE"])
o <- data.frame(indiv=otuN,OTU=rownames(tcount))
o <- sort_df(o,vars="indiv")
ggplot(data=o,aes(x=OTU,y=indiv))+geom_point()+geom_hline(aes(yintercept=0.1*ncol(tcount)),col="red",linetype="dashed")
#+geom_vline(aes(xintercept=o[o$OTU==0.1*ncol(tcount)]),col="red",line="dashed")

#discarted otus analysis
dis <- apply(ones,1,function(x) table(x)["TRUE"]<0.1*ncol(tcount))
dis_df <- tcount[dis,]
ggplot(data=dis_df_m,aes(x=variable,y=value))+geom_point()
dis_df_mean <- apply


#for (i in seq(1,ncol(tcount))){
#  otuN <- apply(ones,1,function(x) table(x)["TRUE"]==i)
#  opi[i,] <- c(i,table(otuN)["TRUE"])
#  i
#}
opi <- opi[seq(1,78),]
ggplot(data=opi,aes(x=indiv,y=OTUs))+geom_point()+scale_y_continuous(limits=c(0,200))+stat_smooth()

###
ones <- tcount[1,]
ones <- apply(tcount,2,function(x) x!=0)
head(ones[,1:5])
o <- apply(ones,1,function(x) table(x)["TRUE"]==1)
head(ones)
table(o)

tcont_o <- tcount[ones,]

ones2 <- rowSums(tcont_o)>ncol(tcont_o)
table(ones2)
tcont_o <- tcont_o[ones2,]

tc <- melt(tcont_o)
summary(tc)
ggplot(data=tc,aes(x=variable,y=value))+geom_point()+stat_smooth()
ggplot(data=tcont_o, aes(x=rownames(tcont_o),y=tcont_o$'117100FOC_G1'))+geom_point()+stat_smooth()


head(ones)
table(ones)
ggplot(data=tcont_o, aes(x=rownames(tcont_o),y=tcont_o$'117100FOC_G1'))+geom_point()

ggplot(data=tcont_o, aes(x=rownames(tcount), fill=..count..))+geom_histogram()





count <- counts[,4:length(counts)]
head(count[,1:5])
cnt <- melt(count)
head(cnt)





#####
k <- sobs.df[,-1]
k[is.na(k)] <- 0
otus <- apply(k,2,max)
summary(otus)
boxplot(otus)
males = as.numeric()
females = as.numeric()

for (i in seq(1:length(otus))){
  #names(otus[1])
  if (grepl("female",names(otus[i]),fixed=T)){
    females <- rbind(females,otus[i])
  }else{
    males <- rbind(males,otus[i])
  }
}
length(males) = length(otus)
length(females) = length(otus)
OTU <- data.frame(all=otus,males=males,females=females)
boxplot(OTU,ylab='No. of OTUs',col=c('red','cadetblue4','darkolivegreen'))
taxa <- taxonomy$Size[taxonomy$Size>1000]
plot(taxa,pch=16,ylab='No. of Reads',xlab='OTU')
sum(taxonomy$Size)

## keep ##

keep <- otu_ab$taxlevel %in% c("2")
y <- subset(otu_ab[1:5],otu_ab$taxlevel %in% c("2"))
y1 <- as.data.frame(y[,c(3,5)])
write.table(y1,file="Phyllums_Taxa.txt",sep="\t")
bp <- ggplot(y1,aes(x='',y=total,fill=taxon,order=-as.numeric(total)))+
  geom_bar(width=1,stat="identity")

pie <- bp+coord_polar("y",start=0)
pie

head(taxonomy)
unknow_reads <- subset(taxonomy,grepl("Bacteria(100);unclassified(100)", taxonomy$Taxonomy,fixed=T))
head(unknow_reads)
sum(unknow_reads$Size)
sum(taxonomy$Size)

# End #
#######



