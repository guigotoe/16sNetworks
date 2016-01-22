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
setwd("/home/torres/Documents/Projects/Metagenome/resutls/all/Networks/")
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
wdeg <- read.table('wdegree.txt',header=T,sep="\t")
ass <- read.table('assortcoef.txt',header=T,sep="\t")
acc <- read.table('avclustcoef.txt',header=T,sep="\t")
asp <- read.table('avshortest.txt',header=T,sep="\t")
bwc <- read.table('betwcent.txt',header=T,sep="\t")
clos <- read.table('closeness.txt',header=T,sep="\t")
eigc <- read.table('eigvectc.txt',header=T,sep="\t")




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
df <- genus_female
title <- "Females"
genus_graph <- function(df,title){
  gm <- as.matrix(prop.table(as.matrix(df[,c(2:(ncol(df)-1))]),1))
  gm <- as.data.frame(gm)
  gm$"id" <- df$"id"
  gm$"age" <- df$"V3"
  agex <- c()
  k=0
  for (i in seq(1:NROW(gf))){  
    if (i==1){
      agex <- append(agex,gf[,ncol(gm)][i])
    }else{
      if(gm[,ncol(gm)][i]== gm[,ncol(gm)][i-1]){
        if (gm[,ncol(gm)][i] %in% agex){k <- k+0.01}
              agex <- append(agex,(gm[,ncol(gm)][i]+k))
      }else{agex <- append(agex,gm[,ncol(gm)][i])} 
    }
  }
  gm$"agex" <- agex
  #colnames(gm)[20]<- "RC9_gut_group"
  colnames(gm)[which(names(gm) == "gut")] <- "RC9_gut_group"
  gm_m <- melt(gm,id.vars=c("id","age","agex"))
  ranks <- unlist(lapply(gm_m$"age", function(x) getrank(x)))
  taxa <- unlist(apply(gm_m,1, function(x) if(as.numeric(x[5])>0.065){(x[4])}else{"xOthers"}))
  gm_m$"ranks" <- as.factor(ranks)
  gm_m$"taxa" <- as.factor(taxa) 
  gmx <- within(gm_m,position <- factor(age,levels=names(sort(table(age)))))
  gmx2 <- aggregate(value ~ id+taxa+age+agex+position+ranks,data=gmx,FUN=sum)
  colourCount = length(unique(gmx2$"taxa"))
  #colourCount=25
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  ggplot(gmx2,aes(x=as.factor(gmx2$agex),y=value,fill=taxa))+scale_shape_discrete(name  ="Genus")+
    geom_bar(with=1,stat="identity")+
    scale_fill_manual(name="Genus",values=colorRampPalette(brewer.pal(8, "Dark2"))(colourCount))+
    scale_x_discrete("Age")+ylab("Proportion")+
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
#* diversity                   *#
#################################
calc_o <- sort_df(calc,vars="nseqs")
for (i in seq(1:NROW(calc_o))){ 
  sample.name <- gregexpr("[[:digit:]]+[[:upper:]]{3}",calc_o$group[i])
  sample.name = regmatches(calc_o$group[i],sample.name)
  calc_o$indiv[i] <- sample.name[[1]][1]
}
head(calc_o)
ggplot(calc_o,aes(x=V3,y=nseqs,color=V3))+geom_point()
sum <- merge(as.data.frame(calc_o),as.data.frame(design),by.x="indiv",by.y="V1")
sum <- sort_df(sum,vars="V3")
ggplot(sum,aes(x=V3,y=nseqs,color=V2))+geom_point()+
  geom_hline(aes(yintercept=mean(sum$nseqs)),col="red",linetype="dashed")+
  xlab("Individuals age")+
  scale_colour_hue(name="Gender",breaks=c("female", "male"),labels=c("Female", "Male"),l=50)  # Legend label, use darker colors
#ggplot(sum,aes(x=V3,y=sobs,col=V2))+geom_point()+ylab("OTUs observed")+xlab("Individuals age")


#Normalized
for (i in seq(1:NROW(calcN))){ 
  sample.name <- gregexpr("[[:digit:]]+[[:upper:]]{3}",calcN$group[i])
  sample.name = regmatches(calcN$group[i],sample.name)
  calcN$indiv[i] <- sample.name[[1]][1]
}
head(calcN)
N <- merge(as.data.frame(calcN),as.data.frame(design),by.x="indiv",by.y="V1")
head(N)
N <- sort_df(N,vars="sobs")
plot(N$sobs)
ggplot(N,aes(x=V3,y=sobs,col=V2))+geom_point()+ylab("OTUs observed")+xlab("Individuals age")+
  scale_colour_hue(name="Gender",breaks=c("female", "male"),labels=c("Female", "Male"),l=50)+  # Legend label, use darker colors
  stat_smooth(method="lm")+ facet_grid(.~V2,margins=TRUE)
ggplot(N,aes(x=V3,y=sobs,col=V2))+geom_point()+ylab("OTUs observed")+xlab("Individuals age")+
  scale_colour_hue(name="Gender",breaks=c("female", "male"),labels=c("Female", "Male"),l=50)+  # Legend label, use darker colors
  stat_smooth(method="loess",formula=y~ns(x,4))+ facet_grid(.~V2,margins=TRUE)

                   

##Alpha diversity normalized
# invsimpson
N$logInvS <- log2(N$invsimpson)
ggplot(N,aes(x=V3,y=invsimpson,col=V2))+geom_point()+ylab("Relative alpha diversity (invSimpson)")+
  xlab("Individuals age")+scale_colour_hue(name="Gender",breaks=c("female", "male"),labels=c("Female", "Male"),l=50)+
  stat_smooth(method="lm")+facet_grid(.~V2,margins=TRUE)



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


#################################
#* Reads in OTU - Distribution *#
#################################

head(counts[,1:5])
tcount <- t(counts[,4:ncol(counts)])
colnames(tcount) <- counts[,2]
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

##################
#* Calculations *#
##################

#############################
## Colector && Rarefaction ##

sobs.df <- OTUs(obs_files) # Colector curve using Observed OTUs

rarefaction$X <- NULL
rare.df <- rarefaction[1]  # Rarefaction curve calculated using Observed OTUs
for (i in seq(from=1,to=length(rarefaction[,-1]),by=3)){rare.df <- cbind(rare.df,rarefaction[,-1][i])}
rare.graph(rare.df)

##Rarefaction normalized
rarefactionN$X <- NULL
rare.df <- rarefactionN[1]  # Rarefaction curve calculated using Observed OTUs
for (i in seq(from=1,to=length(rarefactionN[,-1]),by=3)){rare.df <- cbind(rare.df,rarefactionN[,-1][i])}
#rare.graph(rare.df)
rare.df[,-1] <- log10(rare.df[,-1])
rm <- apply(rare.df[,-1],1,mean)
rstd <- apply(rare.df[,-1],1,sd)
rare.df$mean <- rm
rare.df$std <- rstd
pd <- position_dodge(0.1)
ggplot(rare.df,aes(x=numsampled,y=mean))+
  geom_errorbar(aes(ymin=mean-std,ymax=mean+std),width=.2,size=.3,colour="blue",position=pd)+geom_line(position=pd)+
  geom_point(position=pd, size=2, shape=21, fill="white")+
  xlab("Reads sampled")+ylab("Log10(Observed OTUs)")+ggtitle("Rarefaction or accumulation curve (Normalized sample 11450)")

# End #
#######



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

## calculators ##
boxplot(calc$coverage,col='grey',ylab='Coverage')
boxplot(calc$invsimpson,col='grey',ylab='inv-Simpson')
boxplot(calc$bergerparker,calc$logseries,calc$geometric,calc$bstick,
        calc$bstick,xlab=c('bergerparker','logseries','geometric','bstick'))
boxplot(calc$nseqs,col='grey',ylab='No. of Sequences')

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
#head(rfdata[,1:5])
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
rfdata <- rfdata[noung,]
data <- t(rfdata)
data <- as.data.frame(data)
dim(data)
sum(data)/sum(rdata)
head(data[1:5])
cnphylla <- unlist(lapply(colnames(data),function(x) taxa(x,2)))
cngenus <- unlist(lapply(colnames(data),function(x) taxa(x,6)))

rgenus <- data
colnames(rgenus) <- cngenus
cng <- unique(cngenus)
genus <- data.frame(matrix(NA,nrow=NROW(rgenus),ncol=0))
for (i in 1:length(cng)){
  genus[cng[i]] <- rowSums(rgenus[colnames(rgenus)%in%(cng[i])])
}
genus$"id" <- rownames(rgenus)
genus <- merge(genus,as.data.frame(design),by.x="id",by.y="V1")
genus <- sort_df(genus,vars="V3")

## Males!
genus_male <- subset(genus[,c(1:59,61)],genus$V2=="male")
genus_graph(genus_male,"Males")


gm <- as.matrix(prop.table(as.matrix(genus_male[,c(2:(ncol(genus_male)-1))]),1))
gm <- as.data.frame(gm)
gm$"id" <- genus_male$"id"
gm$"age" <- genus_male$"V3"
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
#colnames(gm)[20]<- "RC9_gut_group"
colnames(gm)[which(names(gm) == "gut")] <- "RC9_gut_group"
gm_m <- melt(gm,id.vars=c("id","age","agex"))
head(gm_m)
ranks <- unlist(lapply(gm_m$"age", function(x) getrank(x)))
taxa <- unlist(apply(gm_m,1, function(x) if(as.numeric(x[5])>0.065){(x[4])}else{"xOthers"}))
head(gm_m)

gm_m$"ranks" <- as.factor(ranks)
gm_m$"taxa" <- as.factor(taxa)

gmx <- within(gm_m,position <- factor(age,levels=names(sort(table(age)))))
head(gmx)
gmx2 <- aggregate(value ~ id+taxa+age+agex+position+ranks,data=gmx,FUN=sum)
head(gmx2)


colourCount = length(unique(gmx2$"taxa"))
#colourCount=25
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

ggplot(gmx2,aes(x=as.factor(gmx2$agex),y=value,fill=taxa))+scale_shape_discrete(name  ="Genus")+
  geom_bar(with=1,stat="identity")+
  scale_fill_manual(name="Genus",values=colorRampPalette(brewer.pal(8, "Dark2"))(colourCount))+
  scale_x_discrete("Age")+ylab("Proportion")+
  guides(fill=guide_legend(ncol=10,keywidth = 0.5, keyheight = 0.5))+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=7),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        legend.position="bottom",legend.box="horizontal",
        legend.text = element_text(size = 4.4),
        legend.title = element_text(size=8, face="bold"))#+facet_wrap(~ ranks,nrow=3,scales="free")#+coord_flip()#

#geom_point(aes(size=value))+
#facet_grid(facets=. ~ ranks)

## Females!

genus_female <- subset(genus[,c(1:59,61)],genus$V2=="female")


genus_graph(genus_female,"Females")


##colnames(sobs.df)[which(names(sobs.df) == "sample.name")] <- sample.name
#length(x) = length(sobs.df$numsampled)
#sobs.df <- cbind(sobs.df,x)