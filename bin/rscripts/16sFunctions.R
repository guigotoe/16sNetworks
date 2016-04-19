####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: 18.04.2016
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

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

##

df <- read.table("/home/torres/Documents/Projects/Metagenome/results/fromMothur/03_16/group_counts",header=F,as.is=T)
gcountg <- function(df){
  plot(sort(df$V3))
  summary(df$V3)
  df[df$V3<13000,]
  
  
}
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

rare.graph <- function(df,lg=TRUE,by.g=FALSE,errbar=FALSE,savef=NULL,title="rarefiPlot"){
  if (lg==TRUE) {
    df[,-1] <- log10(df[,-1])
    lgsc <- ' (log10)'
  }else{lgsc <- ''}
  samples <- melt(df,id.vars="numsampled",na.rm=T)
  if (by.g==TRUE){
    age.group <- apply(samples[2],1, function(x){
      sample.group <- gregexpr("G[[:digit:]]{1}",x) #i
      sample.group = regmatches(x,sample.group) #i
      sample.group = sample.group[[1]][1]
      return(sample.group)
    })
    samples$group <- age.group
    rare.g <- summarySE(samples,measurevar="value",groupvars=c("numsampled","group"),na.rm=FALSE)
    pd <- position_dodge(0.1)
    plot <- ggplot(rare.g,aes(x=numsampled,y=value,colour=group))+geom_line(position=pd)+geom_point(position=pd, size=1, shape=19)+
      xlab("Reads sampled")+ylab(paste("Observed OTUs",lgsc,sep=''))+ggtitle("Rarefaction or accumulation curve")+
      scale_colour_discrete(name = "Age Groups",labels=c("G1(<40)","G2(40,60)","G3(60,80)","G4(80+)"))
    if (errbar == TRUE){
      plot <- plot + geom_errorbar(aes(ymin=value-sd,ymax=value+sd),width=.2,size=.3,colour="blue",position=pd)
      if(!is.null(savef)) ggsave(paste(savef,title,'_errbar.pdf',sep=""),plot=plot,width=12, height=8)
    }else if(!is.null(savef)) ggsave(paste(savef,title,'.pdf',sep=""),plot=plot,width=12, height=8)
      
    }else if (by.g==FALSE){
      plot <- ggplot(samples,aes(x=numsampled,y=value,colour=variable))+geom_line()+xlab("Reads sampled")+
        ylab(paste("Observed OTUs",lgsc,sep=''))+theme(legend.position="none")+ggtitle("Rarefaction or accumulation curve")
      if(!is.null(savef)) ggsave(paste(savef,title,'_by.g.pdf',sep=""),plot=plot,width=12, height=8)
    }
}

##Alpha diversity normalized
# invsimpson and Shannon entropy
div.graphs <- function(df,savef=NULL,title="diversity"){
  df$logInvS <- log2(df$invsimpson)
  df$D1 <- exp(df$shannon)
  
  invsimp <- summarySE(df[,c("age.mean","Gender","logInvS")],measurevar="logInvS",groupvars=c("age.mean","Gender"),na.rm=FALSE)
  ggplot(df[!is.na(df$Gender),],aes(x=Age,y=logInvS,col=Gender))+geom_point()+ylab("Alpha diversity (Log. invSimpson)")+
    xlab("Individuals age")+scale_colour_hue(name="Gender",breaks=c("female", "male"),labels=c("Female", "Male"),l=50)+
    stat_smooth(method="lm")+facet_grid(.~Gender,margins=TRUE)
  if(!is.null(savef)) ggsave(paste(savef,title,'_invSimpson.pdf',sep=""),width=12, height=8)
  ggplot(invsimp,aes(x=age.mean,y=logInvS,col=Gender))+geom_point()+ylab("Alpha diversity (Log. invSimpson)")+
    xlab("Average group's age")+
    scale_colour_hue(name="Gender",l=50)+stat_smooth(method="lm")#+
  #geom_errorbar(aes(ymin=logInvS-sd,ymax=logInvS+sd),width=.2,size=.3,colour="blue")
  if(!is.null(savef)) ggsave(paste(savef,title,'_invSimpsonByG.pdf',sep=""),width=12, height=8)
  D1 <- summarySE(df[,c("age.mean","Gender","D1")],measurevar="D1",groupvars=c("age.mean","Gender"),na.rm=FALSE)
  ggplot(df[!is.na(df$Gender),],aes(x=Age,y=D1,col=Gender))+geom_point()+ylab("Alpha diversity (Shannon entropy)")+
    xlab("Individuals age")+scale_colour_hue(name="Gender",breaks=c("female", "male"),labels=c("Female", "Male"),l=50)+
    stat_smooth(method="lm")+facet_grid(.~Gender,margins=TRUE)
  if(!is.null(savef)) ggsave(paste(savef,title,'_shannonD1.pdf',sep=""),width=12, height=8)
  ggplot(D1,aes(x=age.mean,y=D1,col=Gender))+geom_point()+ylab("Alpha diversity (Shannon entropy)")+
    xlab("Average group's age")+
    scale_colour_hue(name="Gender",l=50)+stat_smooth(method="lm")#+
  #geom_errorbar(aes(ymin=logInvS-sd,ymax=logInvS+sd),width=.2,size=.3,colour="blue")  
  if(!is.null(savef)) ggsave(paste(savef,title,'_shannonD1ByG.pdf',sep=""),width=12, height=8)
  #numseqs
  #ggplot(calc_o,aes(x=V3,y=nseqs,color=V3))+geom_point()
  #sum <- merge(as.data.frame(calc_o),as.data.frame(design),by.x="indiv",by.y="V1")
  #sum <- sort_df(sum,vars="V3")
  #ggplot(sum,aes(x=V3,y=nseqs,color=V2))+geom_point()+
  #  geom_hline(aes(yintercept=mean(sum$nseqs)),col="red",linetype="dashed")+
  #  xlab("Individuals age")+
  #  scale_colour_hue(name="Gender",breaks=c("female", "male"),labels=c("Female", "Male"),l=50)  # Legend label, use darker colors
  #ggplot(sum,aes(x=V3,y=sobs,col=V2))+geom_point()+ylab("OTUs observed")+xlab("Individuals age")
  
  ## calculators ##
  #boxplot(calc$coverage,col='grey',ylab='Coverage')
  #boxplot(calc$invsimpson,col='grey',ylab='inv-Simpson')
  #boxplot(calc$bergerparker,calc$logseries,calc$geometric,calc$bstick,
  #        calc$bstick,xlab=c('bergerparker','logseries','geometric','bstick'))
  #boxplot(calc$nseqs,col='grey',ylab='No. of Sequences')
}



add.GenderAge <- function(df,design,colname="group"){
  df$ID <- rep(NA,NROW(df))
  df$Gender <- rep(NA,NROW(df))
  df$Age <- rep(NA,NROW(df))
  df$age.group <- rep(NA,NROW(df))
  df$age.mean <- rep(NA,NROW(df))
  for (i in seq(1:NROW(df))){
    x <- df[i,]
    y <- design[strsplit(as.character(x[,as.character(colname)]),'_')[[1]][1],]
    df[,"Gender"][df[,as.character(colname)]==x[,as.character(colname)]] <- as.character(y[[1]])
    df[,"Age"][df[,as.character(colname)]==x[,as.character(colname)]] <- as.numeric(y[[2]])
    ID <- gregexpr("[[:digit:]]+[[:upper:]]{3}",as.character(x[,as.character(colname)]))
    ID = regmatches(as.character(x[,as.character(colname)]),ID)
    ID = ID[[1]][1]
    df[,"ID"][df[,as.character(colname)]==x[,as.character(colname)]] <- ID
    if (as.numeric(y[[2]]) <= 40) sample.group <- "G1"
    else if (as.numeric(y[[2]]) > 40 & as.numeric(y[[2]]) <= 60) sample.group <- "G2"
    else if (as.numeric(y[[2]]) > 60 & as.numeric(y[[2]]) <= 80) sample.group <- "G3"
    else if (as.numeric(y[[2]]) > 80) sample.group <- "G4"
    df[,"age.group"][df[,as.character(colname)]==x[,as.character(colname)]] <- sample.group
  }
  for (i in seq(1:NROW(df))){
    x <- df[i,]
    for (j in levels(as.factor(df[,"age.group"]))) {
      if (x[,"age.group"]==j) df[,"age.mean"][df[,as.character(colname)]==x[,as.character(colname)]] <- as.integer(mean(df$"Age"[df$"age.group"==j]))
    }
  }
  return(df)
}
###
# the following give to you the
# minimum OTU bin length. Below that assumed spourious
###
otu.filtering <- function(taxonomy,qval,nval=100,plots=F){
  otu_size_dist <- data.frame("Size"=character(n),"OTUs"= numeric(n),"more_OTUs"= numeric(n),"less_OTUs"= numeric(n),stringsAsFactors = FALSE)
  for (i in seq(1:n)){
    otu_size_dist[i,] <- c(i,table(taxonomy[2]==i)["TRUE"],table(taxonomy[2]>=i)["TRUE"],table(taxonomy[2]<=i)["TRUE"])
  }
  threshold <- quantile(idxs$ace,probs=qval)  ## Estimated population size based on rare OTUs calculated by ace
  otu_size_dist[nval,]$Size <- paste(nval,'+/-',sep=" ")
  min_read_len <- tail(which(otu_size_dist$more_OTUs >= threshold),1)
  #threshold
  #otu_size_dist[which(otu_size_dist$more_OTUs >= threshold),]
  trim_otus <- otu_size_dist[min_read_len:NROW(otu_size_dist),3]
  retained.info <- sum(taxonomy$Size[taxonomy$Size>=as.integer(otu_size_dist$Size[min_read_len])])/sum(taxonomy$Size) # retained info 0.9824082
  if (plots == T){
    require(fitdistrplus)
    
    plot(otu_size_dist[,3],main="Distribution of OTU sizes",ylab="Num. OTUs with more than 'x' Num. sequences",xlab="Num. sequences")
    abline(v=otu_size_dist$Size[min_read_len],col="red",lty=2)
    ##
    plot(trim_otus,main="Distribution of OTU sizes - trimmed data",ylab="Num. OTUs with more than 'x' Num. sequences",xlab="Num. sequences")
    ##
    otu_dist <- fitdist(otu_size_dist[,3],"lnorm")
    plotdist(otu_size_dist[,3],"lnorm",para=list(meanlog=otu_dist$estimate[1],sdlog=otu_dist$estimate[2]))
    otu_trim_dist <- fitdist(trim_otus,"lnorm")
    plotdist(trim_otus,"lnorm",para=list(meanlog=otu_trim_dist$estimate[1],sdlog=otu_trim_dist$estimate[2]))
  }
  return(as.integer(otu_size_dist$Size[min_read_len]))
}
getTaxGenderDF <- function(df.raw,taxlevel,design){
  unique.df <- unlist(unique(df.raw[taxlevel,]))
  df.r <- df.raw[-c((NROW(df.raw)-2):NROW(df.raw)),]
  colnames(df.r) <- df.raw[taxlevel,]
  ## changing all values as.character to numeric -->
  idx <- sapply(df.r,is.character)
  df.r[idx] <- lapply(df.r[idx],function(x) as.numeric(as.character(x)))
  ##
  df <- data.frame(matrix(NA,nrow=NROW(df.r),ncol=0))
  # get sum of all OTUs of the same taxlevel
  for (i in 1:length(unique.df)) df[unique.df[i]] <- rowSums(df.r[colnames(df.r)%in%(unique.df[i])])
  df$"id" <- rownames(df.r)
  rownames(df) <- rownames(df.r)
  df <- add.GenderAge(df,design,colname="id")  # retrieve design info, age, gender, etc.
  df_male <- subset(df,df$Gender=="male") ## Males!
  df_female <- subset(df,df$Gender=="female")  ## Males!
  return(list(all=df,males=df_male,females=df_female))
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

fillzeros <- function(df,output="prop",method="SQ"){
  #output=c("prop","counts"),method=c("GBM","SQ","BL","CZM")
  nodup <- function(x){
    if (is.na(data.frame(table(duplicated(as.factor(x))))$Freq[2] > 0)) return(x[length(x)])
    else nodup(append(x[1:length(x)-1],x[length(x)]+0.01))
  }
  #gm <- as.matrix(prop.table(as.matrix(df[,-c((ncol(df)-5):ncol(df))]),1))  # fractions matrix
  require(zCompositions)
  gm <- df[,-c((ncol(df)-5):ncol(df))]                                      # Compositional Data - proportions
  gm <- cmultRepl(gm,method=method,output=output)
  gm <- as.data.frame(gm)
  gm$"id" <- df$"ID"
  gm$"gender" <- df$"Gender"
  gm$"age" <- df$"Age"
  gm$"age.group" <- df$"age.group"
  gm$"age.mean" <- df$"age.mean"
  agex <- c()
  for (i in gm$age){
    agex <- append(agex,i)
    if (!is.na(data.frame(table(duplicated(as.factor(agex))))$Freq[2] > 0)) agex <- append(agex[1:length(agex)-1],nodup(agex))
  }
  gm$"agex" <- as.numeric(agex)
  return(gm)
}

#df <- phylla.dfs$all
#title <- "Males"
#limit=0.06
#taxlevel="Phylla"
#method="SQ"
#savef=NULL
#ids=F
#output="prop"

tax_graph <- function(df,limit=0.06,taxlevel="Genus",savef=NULL,ids=F,method="SQ"){
  gm <- fillzeros(df,output="prop",method=method)
  colnames(gm)[which(names(gm) == "gut")] <- "RC9_gut_group"
  colnames(gm)[which(names(gm) == "Incertae")] <- "Incertae_Sedis"
  gm_m <- melt(gm,id.vars=c("id","gender","age","agex","age.group","age.mean"))
  ranks <- unlist(lapply(gm_m$"age", function(x) getrank(x)))
  taxa <- unlist(apply(gm_m,1, function(x) if(as.numeric(x["value"])>limit){(x["variable"])}else{"Others"}))
  gm_m$"ranks" <- as.factor(ranks)
  gm_m$"taxa" <- as.factor(taxa) 
  gmx <- within(gm_m,position <- factor(age,levels=names(sort(table(age)))))
  gmx2 <- aggregate(value ~ id+taxa+age+agex+position+ranks+age.group+age.mean+gender,data=gmx,FUN=sum) # Because some taxa now are xOthers so we need to summ their values
  
  ## before plot we need to order the taxa levels according their values of abundance and prevalence in individuals
  taxorder <- data.frame(taxa=levels(gmx2$taxa),rate=rep(0,length(levels(gmx2$taxa))))
  for (i in levels(gmx2$taxa)){
    a <- sum(unlist(lapply(gmx2$value[gmx2$taxa==i],sum))) # OTU global abundace. max = No. of individuals 
    b <- table(gmx2$taxa)[i][[1]]                          # OTU prevalence. max = No. of individuals
    taxorder$rate[taxorder$taxa==i] <- sum(a,b)
  }
  taxorder <- taxorder[order(-taxorder[,"rate"]),]
  taxorder <- taxorder[-which(taxorder$taxa=="Others"),]
  gmx2$taxa <- factor(gmx2$taxa,levels=c(as.character(taxorder$taxa),"Others"))
  gmx2 <- gmx2[order(gmx2$agex,gmx2$taxa),]

  ## placing the colors to the taxas
  colourCount = length(levels(gmx2$taxa))
  if (colourCount <= 10) {palette <- c(colorRampPalette(brewer.pal(8, "Set1"))(colourCount-1),"#3D3D3D")
  }else if (colourCount <= 19) {palette <- c(colorRampPalette(brewer.pal(8, "Set1"))(9),
                                             colorRampPalette(brewer.pal(8, "Accent"))(colourCount-10),"#3D3D3D")
  }else palette <- c(colorRampPalette(brewer.pal(8, "Set1"))(9),colorRampPalette(brewer.pal(8, "Set3"))(9),
                     colorRampPalette(brewer.pal(8, "Accent"))(colourCount-19),"#3D3D3D")
  
  ## Ploting by genders
  # labels for x axis
  
  for (g in levels(as.factor(gm$gender))){
    title <- paste(toupper(substring(as.character(g),1,1)),substring(as.character(g),2),sep="")
    gm.g <- subset(gm,gm$gender==g) ## each gender
    x_labels <- gm.g[,c("agex","age","id")][order(gm.g[,c("agex","age","id")][,"agex"]),]
    if (ids==T){ x_lab <- unlist(apply(x_labels,1,function(x) paste(x[3],x[2],sep="_")))
    }else  x_lab <- x_labels$age
    gmx2.g <- subset(gmx2,gmx2$gender==g) ## each gender
  
    # to be color consistente 
    e <- length(table(gmx2.g$taxa)[table(gmx2.g$taxa)==0])
    if (e!=0) palette.g <- palette[-c((length(palette)-e):(length(palette)-1))]
    else palette.g <- palette
    
    ggplot(gmx2.g, aes(x=as.factor(gmx2.g$agex),y=value,fill=taxa))+
      geom_bar(stat="identity",colour="black")+
      scale_fill_manual(name=taxlevel,values=palette.g) + #limits=c(levels(gmx2[1]),levels(gmx2[length(levels(gmx2))]))
      scale_x_discrete("Age",labels=as.character(x_lab))+ylab("Proportion")+
      guides(fill=guide_legend(ncol=6,keywidth=1, keyheight=1))+
      ggtitle(title)+#facet_wrap(~ age+facet_wrap(~ taxa))+
      theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8),
            panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
            legend.position="bottom",legend.box="horizontal",
            legend.text = element_text(size=10),
            legend.title = element_text(size=12, face="bold"),
            plot.title = element_text(lineheight=.8, face="bold"))
    if(!is.null(savef)) ggsave(paste(savef,title,'_',taxlevel,'.pdf',sep=""),width=12, height=8)
  }
}

change.name <- function(df){
  for (i in seq(1:NROW(df))){ 
    sample.name <- gregexpr("[[:digit:]]+[[:upper:]]{3}",calc_o$group[i])
    sample.name = regmatches(calc_o$group[i],sample.name)
    calc_o$indiv[i] <- sample.name[[1]][1]
  }
}

################
#df.list <- family.dfs
#taxlevel <- "Family"

outliers <- function(df.list,method="SQ",gy.gender=c("male","female")){
  cp <- fillzeros(df.list$all,output="count",method="SQ")
  avec <- c((ncol(cp)-5):ncol(cp))
  cp.x <- cp
  for (i in colnames(cp[,-avec])){
    k <- unlist(lapply(cp[,i],function(x) if(x>2)return(1)))
    if (sum(k)/NROW(cp)<0.1) cp.x[,i] <- NULL
  }
  
  g <- gy.gender
  cp.xg <- subset(cp.x,cp.x$gender%in%g)
  avec.x <- c((ncol(cp.xg)-5):ncol(cp.xg))
  
  #########################
  #   *  Outliers  *   ###
  # robust Mahalanobis distances 
  # internally applies a iso-metric log-ratio transformation to the compositions
  # to search for outliers in the real space.
  
  require(robCompositions)
  
  outliers <- outCoDa(cp.xg[,-avec.x],method="robust")
  cp.xgo <- cbind(cp.xg,mahaDist=outliers$mahalDist)
  return(list(outliers=outliers,cp.xgo=cp.xgo))
  
}

cpanalysis <- function(df.list,taxlevel,savef=NULL,title="cpanalysis"){
  require(robCompositions)
  outl <- outliers(df.list)
  outliers <- outl$outliers
  cp.xgo <- outl$cp.xgo
  
  colourCount = length(levels(as.factor(cp.xgo$age.group)))
  col <- colorRampPalette(brewer.pal(8, "Set1"))(colourCount)
  ggplot(cp.xgo,aes(x=id,y=mahaDist,colour=age.group))+
    geom_point(aes(shape=as.factor(cp.xgo$gender)),size=3.5,colour='black')+
    geom_point(aes(shape=as.factor(cp.xgo$gender)),size=2.5)+
    scale_color_brewer(name="Age Groups",palette="Set3")+
    ggtitle("Outliers detection - threshold quantile 0.975")+
    scale_shape(name="Gender")+
    ylab("Robust Mahalanobis distance")+xlab("Samples")+
    geom_hline(yintercept=outliers$limit,linetype="dashed",color="red")
  
  if(!is.null(savef)) ggsave(paste(savef,title,'_',taxlevel,'_outliers.pdf',sep=""),width=12, height=8)
  
  exc <- cp.xgo[cp.xgo$mahaDist>outliers$limit,]#5,]# df just outliers
  
  ## geting the outliers proportion by age group
  outprop <- data.frame(age.group=levels(as.factor(cp.xgo$age.group)))
  outp.g <- data.frame(age.group=levels(as.factor(cp.xgo$age.group)))
  for (i in levels(as.factor(cp.xgo$age.group))){
    outprop$all[outprop$age.group==i] <- round(NROW(exc[exc$age.group==i,])/NROW(cp.xgo[cp.xgo$age.group==i,]),2)
    outprop$Females.a[outprop$age.group==i] <- round(NROW(exc[exc$age.group==i & exc$gender=="female",])/NROW(cp.xgo[cp.xgo$age.group==i,]),2)
    outprop$Males.a[outprop$age.group==i] <- round(NROW(exc[exc$age.group==i & exc$gender=="male",])/NROW(cp.xgo[cp.xgo$age.group==i,]),2)
    outp.g$Females[outprop$age.group==i] <- round(NROW(exc[exc$age.group==i & exc$gender=="female",])/NROW(cp.xgo[cp.xgo$age.group==i & cp.xgo$gender=="female",]),2)
    outp.g$Males[outprop$age.group==i] <- round(NROW(exc[exc$age.group==i & exc$gender=="male",])/NROW(cp.xgo[cp.xgo$age.group==i & cp.xgo$gender=="male",]),2)
  } 
  outp <- merge(outprop,outp.g,by="age.group",all.x=T)
  
  #outprop.m <- melt(outprop,id.vars="age.group")
  outp.m <- melt(outp,id.vars=c("age.group"))
  outp.m$df <- c(rep("All",12),rep("All by gender",8))
  ggplot(outp.m,aes(x=age.group,y=value,colour=variable,group=variable))+
    geom_point(aes(colour=variable, shape=variable, group=variable),size=3)+
    geom_line(aes(colour=variable, group=variable))+
    xlab("Age group")+ylab("Outlier proportion")+
    scale_colour_discrete(name="Ratio",labels=c("all.o/all.g","oF.g/all.g","oM.g/all.g","oF.g/aF.g","oM.g/aM.g"))+
    scale_shape_discrete(name="Ratio",labels=c("all.o/all.g","oF.g/all.g","oM.g/all.g","oF.g/aF.g","oM.g/aM.g"))+
    ggtitle("Proportion of outliers by age group")+facet_wrap(~df)
  
  if(!is.null(savef)) ggsave(paste(savef,title,'_',taxlevel,'_outliersProp.pdf',sep=""),width=12, height=8)
  
  ### confirmation with replicated libraries ####
  dups <- read.table(paste(files.path,'duplicated.txt',sep=''),header=F,sep="\t")
  dups$counter <- rep(0,NROW(dups))
  dups$exc <- rep(NA,NROW(dups))
  for (j in rownames(exc)) {
    i <- strsplit(j,'_')[[1]][1]
    if(i%in%dups[,1]){
      dups$counter[dups$V1==i] <- dups$counter[dups$V1==i]+1
      dups$exc[dups$V1==i] <- j
    }else if(i%in%dups[,2]){
      dups$counter[dups$V2==i] <- dups$counter[dups$V2==i]+1
      dups$exc[dups$V2==i] <- j
    }
  }
  
  dup.exc.prc <- c(Prop.both.libs=round(1-(length(dups$counter[dups$counter==1])/NROW(dups)),2),
                   Num.both.libs=(NROW(dups)-length(dups$counter[dups$counter==1])),Num.libs=NROW(dups))
  print(dup.exc.prc)
  #####
  #cp.yo = cp.xgo[!outliers$outlierIndex,]
  #cp.y = cp.ya
  #cp.ya = cp.xgo
  cp.y <- cp.xgo
  avec.y <- c((ncol(cp.y)-6):ncol(cp.y))
  pca <- pcaCoDa(cp.y[-avec.y],method="robust")
  pca$explainedVar <- pca$princompOutputClr$sdev^2/sum(pca$princompOutputClr$sdev^2)
  ggplot(as.data.frame(pca$scores),aes(x=Comp.1,y=Comp.2,colour=cp.y$age.group))+
    geom_hline(yintercept = 0,colour="gray65")+geom_vline(xintercept = 0,colour="gray65")+
    geom_point(aes(colour=cp.y$age.group, shape=cp.y$gender,group=cp.y$age.group),size=3)+
    scale_colour_brewer(name="Age group",palette="Set3")+
    scale_shape_discrete(name=taxlevel)+
    xlab(paste("PCA1 (clr-robust) ",round(pca$explainedVar[1]*100,2),"%",sep=''))+
    ylab(paste("PCA2 (clr-robust) ",round(pca$explainedVar[2]*100,2),"%",sep=''))+
    ggtitle(paste("PCA plot of idividuals - by ",taxlevel,sep=''))
  
  if(!is.null(savef)) ggsave(paste(savef,title,'_',taxlevel,'_pca.pdf',sep=""),width=12, height=8)
  
  ## Circle of correlations
  
  # create a circle of ratio 1
  circle <- function(center=c(0,0),npoints=100){
    r = 1
    tt = seq(0,2*pi,length=npoints)
    xx = center[1]+r*cos(tt)
    yy = center[1]+r*sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  
  corcir = circle(c(0, 0), npoints = 100) 
  correlations = as.data.frame(cor(cp.y[-avec.y], pca$scores))
  arrows = data.frame(x1 = rep(0,NROW(correlations)), y1 = rep(0,NROW(correlations)), x2 = correlations$Comp.1, 
                      y2 = correlations$Comp.2)
  
  ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
    geom_segment(data=arrows, aes(x = x1, y = y1, xend = x2, yend = y2),colour="black",
                 arrow = arrow(length = unit(0.03, "npc"))) +
    geom_point(data = arrows,aes(x = x2, y = y2))+
    geom_text_repel(data=correlations, aes(x = Comp.1, y = Comp.2, label = rownames(correlations))) + 
    geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0,colour = "gray65") + 
    xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "PCA1 aixs",y = "PCA2 axis") + 
    ggtitle("Circle of correlations")

  if(!is.null(savef)) ggsave(paste(savef,title,'_',taxlevel,'_corrcircle.pdf',sep=""),width=8, height=8)
  ## Discriminant analysis by Fisher rule
  #y <- cp.y[-avec.y]
  #grp<- cp.y$age.group
  #daf <- daFisher(y,grp=grp,method="robust")
}

#########
#* END *#
#########
