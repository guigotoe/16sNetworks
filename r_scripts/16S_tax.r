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
library(reshape)
library(reshape2)
library(scales)
library(vegan)
library(ade4)
library(RColorBrewer) 


#* Globals *#
script_loc = getwd()
setwd("/home/torres/Documents/Projects/Metagenome/resutls/all/graph_tax/")
files.path = '/home/torres/ikmb_storage/Metagenome/16s/all/last/'#
design <- read.table(paste(files.path,'design.txt',sep=''),header=F,sep="\t")
rownames(design) <- design$V1 
obs_files <- list.files(files.path,pattern="\\.sobs$")
rarefaction <- read.table(paste(files.path,'16S.an.groups.rarefaction',sep=''),header=T,sep="\t")
rarefactionN <- read.table(paste(files.path,'16S.an.0.03.subsample.groups.rarefaction',sep=''),header=T,sep="\t")
taxonomy <- read.table(paste(files.path,'16S.an.0.03.cons.taxonomy',sep=''),header=T,sep="\t")
calc <- read.table(paste(files.path,'16S.an.groups.summary',sep=''),header=T,sep="\t")
calcN <- read.table(paste(files.path,'16S.an.0.03.subsample.groups.summary',sep=''),header=T,sep="\t")
otu_ab <- read.table(paste(files.path,'16S.an.cons.taxonomy',sep=''),header=T,sep="\t")
counts <- read.table(paste(files.path,'16S.an.0.03.subsample.shared',sep=''),header=T,sep="\t") 

###############
#* Functions *#
###############

OTUs <- function (file_list){
  df = data.frame(numsampled=1)
  for (i in seq(1:length(file_list))){ 
    x <- read.table(paste(files.path,file_list[i],sep=''),header=T,sep="\t") #i 10-1-100
    sample.name <- gregexpr("[[:digit:]]+[[:upper:]]{3}",file_list[i]) #i
    if (sample.name[[1]][1] != -1){
      sample.name = regmatches(file_list[i],sample.name) #i
      sample.name = sample.name[[1]][1]
      sample.name = paste(sample.name,design[sample.name,2],design[sample.name,3],sep="_")
    }else{sample.name = tail(strsplit(file_list[i],"\\.")[[1]],n=2L)[1]}
    colnames(x)[2] <- sample.name
    df <- merge(df,as.data.frame(x),by='numsampled',all.x=T,all.y=T)
  }
  return(df)
}
rare.graph <- function(df){
  df[,-1] <- log10(df[,-1])
  samples <- melt(df,id.vars="numsampled",na.rm=T)
  ggplot(samples,aes(x=numsampled,y=value,colour=variable))+geom_line()+xlab("Reads sampled")+
    ylab("Log10(Observed OTUs)")+theme(legend.position="none")+ggtitle("Rarefaction or accumulation curve") 
}

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
brewerplot <- function (palette) {
  p + scale_fill_brewer(palette = palette) + opts(title=palette)
}
taxa <- function(tax,n){ 
  tax <- strsplit(as.character(taxonomy$Taxonomy[taxonomy$OTU==tax][1]),";")
  return(gsub("[[:punct:]]","",unlist(regmatches(tax[[1]][n],gregexpr("[[:alpha:]]+[[:punct:]]",tax[[1]][n])))[1]))
}
indiv <- function(x){
  sample.name <- gregexpr("[[:digit:]]+[[:upper:]]{3}",x)
  sample.name = regmatches(x,sample.name)
  return(sample.name[[1]][1])
}
getrank <- function (i){
  if (i %in% c(0:19)){return("< 20")
  }else if (i %in% c(20:29)){return("20-30")
  }else if (i %in% c(30:34)){return("30-35")
  }else if (i %in% c(35:39)){return("35-40")
  }else if (i %in% c(40:44)){return("40-45")
  }else if (i %in% c(45:49)){return("45-50")
  }else if (i %in% c(50:54)){return("50-55")
  }else if (i %in% c(55:59)){return("55-60")
  }else if (i %in% c(60:64)){return("60-65")
  }else if (i %in% c(65:69)){return("65-70")
  }else if (i %in% c(70:74)){return("70-75")
  }else if (i %in% c(75:79)){return("75-80")
  }else if (i %in% c(80:84)){return("80-85")
  }else if (i %in% c(85:90)){return("85-90")
  }else if (i %in% c(91:130)){return("90+")
  }else {return(NA)}
}
df <- genus_male
df <- phyll_male
title <- "Genus - Males"
title <- "Phyllums - Males"
genus_graph <- function(df,title){
  gm <- as.matrix(prop.table(as.matrix(df[,c(2:(ncol(df)-1))]),1))  # get the proportions
  gm <- as.data.frame(gm)
  gm$"id" <- df$"id"
  gm$"age" <- df$"V3"
  agex <- c()
  k=0
  for (i in seq(1:NROW(gm))){  
    if (i==1){
      agex <- append(agex,gm[,ncol(gm)][i])
    }else{
      if(gm[,ncol(gm)][i]== gm[,ncol(gm)][i-1]){
        if (gm[,ncol(gm)][i] %in% agex){k <- k+0.01}
        agex <- append(agex,(gm[,ncol(gm)][i]+k))
      }else{agex <- append(agex,gm[,ncol(gm)][i])} 
    }
  }
  gm$"agex" <- agex
  colnames(gm)[which(names(gm) == "gut")] <- "RC9_gut_group"
  gm_m <- melt(gm,id.vars=c("id","age","agex"))
  ##gm_m <- melt(gm,id.vars=c("id","age"))
  ranks <- unlist(lapply(gm_m$"age", function(x) getrank(x)))
  taxa <- unlist(apply(gm_m,1, function(x) if(as.numeric(x[ncol(gm_m[1,])])>0.065){(x[ncol(gm_m[1,])-1])}else{"xOthers"}))
  gm_m$"ranks" <- as.factor(ranks)
  #gm_m$"taxa" <- as.factor(gm_m$"variable") # Just for phyllums
  gm_m$"taxa" <-as.factor(taxa) # Just for genders
  gmx <- within(gm_m,position <- factor(age,levels=names(sort(table(age)))))
  gmx2 <- aggregate(value ~ id+taxa+age+agex+position+ranks,data=gmx,FUN=sum)
  #gmx2 <- aggregate(value ~ id+taxa+age+position,data=gmx,FUN=sum)
  colourCount = length(unique(gmx2$"taxa"))
  #colourCount=8
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  #as.factor(gmx2$agex)
  ggplot(gmx2,aes(x=as.factor(gmx2$agex),y=value,fill=taxa))+scale_shape_discrete(name  ="Genus")+
    geom_bar(with=1,stat="identity")+
    scale_fill_manual(name="Genus",values=colorRampPalette(brewer.pal(8, "Dark2"))(colourCount))+
    scale_x_discrete("Age",breaks=gmx2$'agex',labels=gmx2$'age')+
    ylab("Proportion")+
    guides(fill=guide_legend(ncol=10,keywidth = 0.5, keyheight = 0.5))+
    ggtitle(title)+
    theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=7),
          panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          legend.position="bottom",legend.box="horizontal",
          legend.text = element_text(size = 6),
          legend.title = element_text(size=8, face="bold"),
          plot.title = element_text(lineheight=.8, face="bold"))#+facet_wrap(~ ranks,nrow=3,scales="free")#+coord_flip()#
}
#########
#* END *#
#########

#################################
#* Taxonomy COA                *#
#################################

##Filtering ##
head(counts[,1:5])
rdata <- counts[,4:(ncol(counts)-1)]
indivs <- unlist(lapply(counts[,2],function(x) indiv(x)))
rownames(rdata) <- indivs
#head(rdata[,(ncol(rdata)-4):ncol(rdata)])
abotus <- apply(rdata,1,function(x) x>=10)
#head(abotus[,1:5])
abotusI <- apply(abotus,1,function(x) table(x)["TRUE"]>=5)
abotusI[is.na(abotusI)] <- FALSE
#table(abotusI)
#head(abotusI)
rfdata <- t(rdata)
rfdata <- rfdata[abotusI,]
head(rfdata[,1:5])
dim(rfdata)
sum(rfdata)/sum(rdata)
phylla_all <- unlist(lapply(rownames(rfdata),function(x) taxa(x,2)))
genus_all <- unlist(lapply(rownames(rfdata),function(x) taxa(x,6)))
#head(phylla)
noun <- phylla_all!="unclassified"
noung <- genus_all!="unclassified"
#head(noun)

##** GENUS **##
## Removing unknown genus
data <- t(rfdata[noung,])  # filtering the unclassified)
data <- as.data.frame(data)
dim(data)
sum(data)/sum(rdata)
head(data[1:5])
cnphylla <- unlist(lapply(colnames(data),function(x) taxa(x,2)))
cngenus <- unlist(lapply(colnames(data),function(x) taxa(x,6)))

rgenus <- data
colnames(rgenus) <- cngenus
head(rgenus[1:5])
cng <- unique(cngenus)
genus <- data.frame(matrix(NA,nrow=NROW(rgenus),ncol=0))
for (i in 1:length(cng)){
  genus[cng[i]] <- rowSums(rgenus[colnames(rgenus)%in%(cng[i])])
}
genus$"id" <- rownames(rgenus)
head(genus[1:5])
genus <- merge(genus,as.data.frame(design),by.x="id",by.y="V1")
genus <- sort_df(genus,vars="V3")

## Males!

genus_male <- subset(genus,genus$V2=="male")
genus_male <- genus_male[,!(colnames(genus_male)=="V2")]
genus_graph(genus_male,"Males")

## Females!

genus_female <- subset(genus,genus$V2=="female")
genus_female <- genus_female[,!(colnames(genus_female)=="V2")]
genus_graph(genus_female,"Females")

##** Phyllum **##

## Removing unknew taxlevel
data <- t(rfdata) # filtering the unclassified
data <- as.data.frame(data)
dim(data)
sum(data)/sum(rdata)
head(data[1:5])
cnphylla <- unlist(lapply(colnames(data),function(x) taxa(x,2)))
#cngenus <- unlist(lapply(colnames(data),function(x) taxa(x,6)))

colnames(data) <- cnphylla
head(data[1:5])
ucol <- unique(cnphylla)
phyll<- data.frame(matrix(NA,nrow=NROW(data),ncol=0))
for (i in 1:length(ucol)){
  phyll[ucol[i]] <- rowSums(data[colnames(data)%in%(ucol[i])])
}
phyll$"id" <- rownames(data)
head(phyll[1:5])
phyll <- merge(phyll,as.data.frame(design),by.x="id",by.y="V1")
phyll <- sort_df(phyll,vars="V3")

## Males!

phyll_male <- subset(phyll,phyll$V2=="male")
phyll_male <- phyll_male[,!(colnames(phyll_male)=="V2")]
genus_graph(genus_male,"Males")

## Females!

phyll_female <- subset(phyll,phyll$V2=="female")
phyll_female <-phyll_female[,!(colnames(phyll_female)=="V2")]
genus_graph(genus_female,"Females")


##colnames(sobs.df)[which(names(sobs.df) == "sample.name")] <- sample.name
#length(x) = length(sobs.df$numsampled)
#sobs.df <- cbind(sobs.df,x)



################# networks ##########
###############################################################
### Analysis by Age groups

df.g1 <- data[meta$age.group=="G1",]
df.g2 <- data[meta$age.group=="G2",]
df.g3 <- data[meta$age.group=="G3",]
df.g4 <- data[meta$age.group=="G4",]
df.f <- data[meta$Gender=="female",]
df.f.g1 <- data[meta$Gender=="female"&meta$age.group=="G1",]
df.f.g2 <- data[meta$Gender=="female"&meta$age.group=="G2",]
df.f.g3 <- data[meta$Gender=="female"&meta$age.group=="G3",]
df.f.g4 <- data[meta$Gender=="female"&meta$age.group=="G4",]
df.m <- data[meta$Gender=="male",]
df.m.g1 <- data[meta$Gender=="male"&meta$age.group=="G1",]
df.m.g2 <- data[meta$Gender=="male"&meta$age.group=="G2",]
df.m.g3 <- data[meta$Gender=="male"&meta$age.group=="G3",]
df.m.g4 <- data[meta$Gender=="male"&meta$age.group=="G4",]

df.g1.mc <- most_common(df.g1,p=60)
df.f.g1.mc <- most_common(df.f.g1,p=60)
df.m.g1.mc <- most_common(df.m.g1,p=60)

se.mb.dfg1mc <- spiec.easi(df.g1.mc$data,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=50))
ig.mb.dfg1mc <- graph.adjacency(se.mb.dfg1mc$refit, mode='undirected')
vsize.dfg1mc <- rowMeans(clr(df.g1.mc$data, 1))+6
am.coord.dfg1mc <- layout.fruchterman.reingold(ig.mb.dfg1mc)

se.mb.dffg1mc <- spiec.easi(df.f.g1.mc$data,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=50))
ig.mb.dffg1mc <- graph.adjacency(se.mb.dffg1mc$refit, mode='undirected')
vsize.dffg1mc <- rowMeans(clr(df.f.g1.mc$data, 1))+6
am.coord.dffg1mc <- layout.fruchterman.reingold(ig.mb.dffg1mc)

se.mb.dfmg1mc <- spiec.easi(df.m.g1.mc$data,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=50))
ig.mb.dfmg1mc <- graph.adjacency(se.mb.dfmg1mc$refit, mode='undirected')
vsize.dfmg1mc <- rowMeans(clr(df.g1.mc$data, 1))+6
am.coord.dfmg1mc <- layout.fruchterman.reingold(ig.mb.dfmg1mc)


par(mfrow=c(1,3))
plot(ig.mb.dfg1mc, layout=am.coord.dfg1mc, vertex.size=vsize.dfg1mc, vertex.label=NA, main="MB-G1")
plot(ig.mb.dffg1mc, layout=am.coord.dffg1mc, vertex.size=vsize.dffg1mc, vertex.label=NA, main="MB-G1.f")
plot(ig.mb.dfmg1mc, layout=am.coord.dfmg1mc, vertex.size=vsize.dfmg1mc, vertex.label=NA, main="MB-G1.m")

plot(ig.mb.dfg1mc, vertex.size=vsize, vertex.label=NA, main="MB")
dd.mb <- degree.distribution(ig.mb.dfg1mc)

plot(0:(length(dd.mb)-1), dd.mb, ylim=c(0,.35), type='b',col="red")#, 
#ylab="Frequency", xlab="Degree", main="Degree Distributions")
legend("topright", c("MB"), col=c( "red"), pch=1, lty=1)


comdata <- most_common(rdata,p=30)
densityplot(comdata$presence)
densityplot(comdata$presence.keept)
dim(comdata$data)

#data("amgut1.filt")
#head(amgut1.filt)
#dim(amgut1.filt)
#depths <- rowSums(amgut1.filt)

comgut <- comdata$data
depths.x <- rowSums(as.matrix(comgut))


se.mb.gergut <- spiec.easi(comgut,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=50))
## Create igraph objects
ig.mb <- graph.adjacency(se.mb.gergut$refit, mode='undirected')
ig <- graph.adjacency(se.mb.gergut$refit, atrr="weight",mode='undirected')
palf <- colorRampPalette(c("gold","darkorange"))
heatmap(ig.mb[,17:1], Rowv = NA, Colv = NA, col = palf(100), 
        scale="none", margins=c(10,10) )

colnames(ig.mb) <- colnames(se.mb.gergut$data)

## set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(comgut, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)
#par(mfrow=c(1,3))
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.mb, vertex.size=vsize, vertex.label=NA, main="MB")
dd.mb <- degree.distribution(ig.mb)

plot(0:(length(dd.mb)-1), dd.mb, ylim=c(0,.35), type='b',col="red")#, 
#ylab="Frequency", xlab="Degree", main="Degree Distributions")
legend("topright", c("MB"), col=c( "red"), pch=1, lty=1)

## threshold based on degree

max(dd.mb)
min(dd.mb)
mean(dd.mb)

