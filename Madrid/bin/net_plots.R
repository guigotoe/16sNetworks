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
library(ade4)
library(RColorBrewer) 
library(plyr)


#* Globals *#
script_loc = getwd()
setwd("~/Documents/Projects/Metagenome/Madrid/results/")
files.path = '~/Documents/Projects/Metagenome/Madrid/bin/'#
prefix <- "JP_"



m_wdeg <- read.table(paste(files.path,'OTU_best5p/OTUp_males_best5p_wdegree_st.txt',sep=''),header=T,sep="\t")
f_wdeg <- read.table(paste(files.path,'OTU_best5p/OTUp_females_best5p_wdegree_st.txt',sep=''),header=T,sep="\t")

m_ass <- read.table(paste(files.path,'males_assortcoef.txt',sep=''),header=T,sep="\t")
f_ass <- read.table(paste(files.path,'females_assortcoef.txt',sep=''),header=T,sep="\t")

m_bwc <- read.table(paste(files.path,'males_betwcent.txt',sep=''),header=T,sep="\t")
f_bwc <- read.table(paste(files.path,'females_betwcent.txt',sep=''),header=T,sep="\t")

m_clos <- read.table(paste(files.path,'males_closeness.txt',sep=''),header=T,sep="\t")
f_clos <- read.table(paste(files.path,'females_closeness.txt',sep=''),header=T,sep="\t")

m_eigc <- read.table(paste(files.path,'males_eigvectc.txt',sep=''),header=T,sep="\t")
f_eigc <- read.table(paste(files.path,'females_eigvectc.txt',sep=''),header=T,sep="\t")

m_acc <- read.table(paste(files.path,'males_avclustcoef.txt',sep=''),header=T,sep="\t")
f_acc <- read.table(paste(files.path,'females_avclustcoef.txt',sep=''),header=T,sep="\t")

asp <- read.table('avshortest.txt',header=T,sep="\t")


##################
##* Functionns *##
##################

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

df <- m_wdeg
sample <- "Males"

OTUstGraph <- function (df,sample){
  z = data.frame("X"=NULL,"N"=NULL,"P"=NULL,"V"=NULL)
  for (i in seq(1,ncol(t(df)))){
    x = t(df[i,])
    f <- sort(x[2:length(x)],decreasing=T)[1]
    fn <- rownames(x)[which(x %in% f)[1]]
    s <- sort(x[2:length(x)],decreasing=T)[2]
    sn <- rownames(x)[which(x %in% s)[1]]
    t <- sort(x[2:length(x)],decreasing=T)[3]
    tn <- rownames(x)[which(x %in% t)[1]]
    l <- sort(x[2:length(x)],decreasing=F)
    l <- l[min( which ( l != 0 )) : max( which( l != 0 ))][1]
    ln <- rownames(x)[match(x[grep(l,x)],x)][1]
    k = data.frame("X"=rep(x[1],4),"N"=c(fn,sn,tn,ln),"P"=c("F","S","T","L"),"V"=c(f,s,t,l))
    z <- rbind(z,k)
  }
  maxS <- sort(z$V[z$P=="F"],decreasing=T)[1]
  maxS <- which(z$V %in% maxS)[1]
  z[maxS,1]
  minS <- sort(z$V[z$P=="F"],decreasing=F)[1]
  maxS <- which(z$V %in% minS)[1]
  z[minS,1]
  
  
  m <- melt(df,id.vars=c("X"))
  pdf(paste("Best5p_TopOTUstrength",sample,".pdf",sep=""), pointsize = 12, width = 10 , height = 7 ) 
  ggplot(data=m)+geom_line(aes(x=z$X,y=z$V,group=z$P,colour=z$P))+geom_point(aes(x=z$X,y=z$V,colour=z$P),size=1.5)+
    ylab("Strength")+xlab("Window age mean")+ggtitle(paste(sample," - OTU with highest strength",sep=""))+#scale_y_discrete(labels="")+
    scale_colour_discrete(name="Positions",breaks=c("F", "S","T", "L"),labels=c("First","Second","Thirth","Last"))
  dev.off()
  
  pdf(paste("Best5p_AllOTUstrength",sample,".pdf",sep=""), pointsize = 12, width = 10 , height = 7 ) 
  ggplot(data=m,aes(x=X,y=value,group=variable,colour=variable))+geom_line()+geom_point(size=1.5)+
    ylab("Strength")+xlab("Window age mean")+theme(legend.position="none")+ggtitle(paste(sample," - Strength by OTU",sep=""))
  dev.off()
}

#################################
#* Networks metrics graphs     *#
#################################

#-----------#
#  Streght  #
#-----------#

head(m_wdeg[,1:5])
qplot(x=m_wdeg$'X',y=m_wdeg[,4],data=m_wdeg,geom="line")
qplot(x=f_wdeg$'X',y=f_wdeg[,4],data=f_wdeg,geom="line")

#m_wdegx <- m_wdeg(na.omit(m_wdeg))
m_wdm <- apply(m_wdeg[2:ncol(m_wdeg)],1,mean)
m_wdegt <- melt(m_wdeg,id.vars=c("X"))
m_wdSE <- summarySE(m_wdegt,measurevar="value",groupvars=c("X"),na.rm=TRUE)
m_wdSE$"gender" <- rep("male",NROW(m_wdSE))

#f_wdeg[is.na(f_wdeg)] <- 0
f_wdm <- apply(f_wdeg[2:ncol(f_wdeg)],1,mean)
f_wdegt <- melt(f_wdeg,id.vars=c("X"))
f_wdSE <- summarySE(f_wdegt,measurevar="value",groupvars=c("X"),na.rm=TRUE)
f_wdSE$"gender" <- rep("female",NROW(f_wdSE))

wdSE <- rbind(m_wdSE,f_wdSE)

pdf("Best5p_strength.pdf", pointsize = 12, width = 10 , height = 7 ) 
ggplot(data=wdSE,aes(X,value,color=gender))+geom_line()+geom_point()+
  ylab("Mean of strength")+xlab("Window age mean")+
  ggtitle("Strength")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1)
dev.off()
colnames(wdSE) <- c("WinAge","NumOTUs","aveStrength","stdDev","stdError","ConfInterv","Gender")
write.table(wdSE,file="Best5p_strength_Stats.txt",sep="\t",row.names=FALSE)

#-------------#
#  Closeness  #
#-------------#

head(m_clos[,1:5])
qplot(x=m_clos$'X',y=m_clos[,4],data=m_clos,geom="line")
qplot(x=f_clos$'X',y=f_clos[,4],data=f_clos,geom="line")

m_clt <- melt(m_clos,id.vars=c("X"))
m_clSE <- summarySE(m_clt,measurevar="value",groupvars=c("X"))
m_clSE$"gender" <- rep("male",NROW(m_clSE))

f_clt <- melt(f_clos,id.vars=c("X"))
f_clSE <- summarySE(f_clt,measurevar="value",groupvars=c("X"))
f_clSE$"gender" <- rep("female",NROW(f_clSE))

clSE <- rbind(m_clSE,f_clSE)

pdf("JP_closeness.pdf", pointsize = 12, width = 10 , height = 7 ) 
ggplot(data=clSE,aes(X,value,color=gender))+geom_line()+geom_point()+
  ylab("Mean of closeness")+xlab("Window age mean")+
  ggtitle("Closeness")+geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1)
dev.off()


#--------------------------#
#  Eigenvector Centrality  #
#--------------------------#

head(m_eigc[,1:5])
qplot(x=m_eigc$'X',y=m_eigc[,4],data=m_eigc,geom="line")
qplot(x=f_eigc$'X',y=f_eigc[,4],data=f_eigc,geom="line")

m_eit <- melt(m_eigc,id.vars=c("X"))
m_eiSE <- summarySE(m_eit,measurevar="value",groupvars=c("X"))
m_eiSE$"gender" <- rep("male",NROW(m_eiSE))

f_eit <- melt(f_eigc,id.vars=c("X"))
f_eiSE <- summarySE(f_eit,measurevar="value",groupvars=c("X"))
f_eiSE$"gender" <- rep("female",NROW(f_eiSE))

eiSE <- rbind(m_eiSE,f_eiSE)

pdf("JP_eigenvector.pdf", pointsize = 12, width = 10 , height = 7 ) 
ggplot(data=eiSE,aes(X,value,color=gender))+geom_line()+geom_point()+
  ylab("Mean of Eigenvector Centrality")+xlab("Window age mean")+
  ggtitle("Eigenvector Centrality")+geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1)
dev.off()

#----------------------------------#
#  Average clustering coefficient  #
#----------------------------------#

head(m_acc)
m_accx <- m_acc
f_accx <- f_acc

m_accx$"gender" <- rep("male",NROW(m_acc))
f_accx$"gender" <- rep("female",NROW(f_acc))
acc <- rbind(m_accx,f_accx)

pdf("JP_assortativity.pdf", pointsize = 12, width = 10 , height = 7 ) 
ggplot(data=acc,aes(X,ave_cluster_coef,color=gender))+geom_line()+geom_point()+
  ylab("Average assortativity")+xlab("Window age mean")+
  ggtitle("Assortativity coefficient")#+geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1)
dev.off()

#------------------------------#
#  Assortativity coefficient  #
#------------------------------#

head(m_ass)
m_assx <- m_ass
f_assx <- f_ass

m_assx$"gender" <- rep("male",NROW(m_ass))
f_assx$"gender" <- rep("female",NROW(f_ass))
ass <- rbind(m_assx,f_assx)

pdf("JP_clustering.pdf", pointsize = 12, width = 10 , height = 7 ) 
ggplot(data=ass,aes(X,ave_cluster_coef,color=gender))+geom_line()+geom_point()+
  ylab("Average clustering coef.")+xlab("Window age mean")+#ylim(-0.004,-0.003)+
  ggtitle("Clustering coefficient")#+geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1)
dev.off()

#----------------------------------#
#  working with OTU_profile (cor)  #
#----------------------------------#
OTUm <- read.table(paste(files.path,'OTU_best5p/OTUp_males_35.2_best5p.txt',sep=''),header=T,sep="\t")
OTUf <- read.table(paste(files.path,'OTU_best5p/OTUp_females_35.4_best5p.txt',sep=''),header=T,sep="\t")

OTU <- OTUf
#head(OTU[,1:5])
m <- upper.tri(as.matrix(OTU[2:ncol(OTU)]))
up <- OTU[2:ncol(OTU)][m]
up_x <- sort(up)
up_x <- up_x[min( which ( up_x != 0 )) : max( which( up_x != 0 ))]
pdf("Best5p_hisCosineDistance_Female_35-4.pdf", pointsize = 12, width = 10 , height = 7 )
hist(up_x,xlab='Cosine Distance',main="Females cosine distance distribution (2242 links) window age = 35.4")
dev.off()


hist(up)
d <- density(up)
plot(d, main="Kernel Density")
polygon(d, col="red", border="blue")
#boxplot(up,pch=20)
max(up)
min(up)


d_x <- density(up_x)
plot(d_x, main="Kernel Density")
polygon(d_x, col="red", border="blue")

m_wdeg <- read.table(paste(files.path,'OTUp_t01/OTUp_males_t01_wdegree.txt',sep=''),header=T,sep="\t")
f_wdeg <- read.table(paste(files.path,'OTUp_test/OTUp_10test_females_t01_wdegree.txt',sep=''),header=T,sep="\t")

OTUstGraph(m_wdeg,"Males")
OTUstGraph(f_wdeg,"Females")


##** exporting nets to be ploted in cytoscape **##
c=1
for(i in c('OTUp_t01/OTUp_males_33.3_t01.txt','OTUp_t01/OTUp_males_57.5_t01.txt')){
  matrix <- read.table(paste(files.path,'OTUp_t01/OTUp_males_33.3_t01.txt',sep=''),header=T,sep="\t")

  net <- melt(matrix,id.vars=c("X"))
  net[net == 0] <- NA
  row.nolink <- apply(net, 1, function(x){any(is.na(x))})
  net.filtred <- net[!row.nolink,]
  
  x <- apply(net.filtred,1,function(x){
    paste(x[1],'-',x[2],sep="")
    })
  y <- apply(net.filtred,1,function(x){
    paste(x[2],'-',x[1],sep="")
  })
  net.filtred$R <- x
  net.filtred$L <- y
  head(net.filtred)
  s=1
  tf <- apply(net.filtred,1,function(x){
    if (x[5] == net.filtred[s,4] & x[4] == net.filtred[s,5]){
      s = s+1
      return(FALSE)
      print s
      }else {return(TRUE)}
  })
  
  #net.x <- subset(net.filtred, !duplicated(L))[,1:3]
  write.table(net.filtred,file=paste("net",c,".txt",sep=""),sep="\t",col.names=F,quote=F,row.names=F)
  c <- c+1
}

