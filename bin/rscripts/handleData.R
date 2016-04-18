file.dir <- '/home/torres/Documents/Projects/Metagenome/docs/16S'
femke <- read.table(paste(file.dir,'first_info_Femke/FoCus_BSPSPC_Age-Gender_femke.txt',sep = '/'),header = T,sep = "\t",na.strings = "NA")
data.2015 <- read.table(paste(file.dir,'design.txt',sep='/'),header=T,sep="\t")
#last49 <- read.table(paste(file.dir,'last49/16S_Last49.txt',sep='/'),header=T,sep="\t")
#samples <- read.table('/home/torres/Documents/Projects/Metagenome/data/16s/samples.txt',header=F,sep="\t")


femke <- femke[complete.cases(femke[,c("age.years","gender")]),]
femke.nonused <- subset(femke, !femke$new_id%in%data.2015$ID)
g1.f <- femke.nonused[sample(rownames(femke.nonused[femke.nonused$gender=="female" & femke.nonused$age.years>19 & femke.nonused$age.years<=40,]), 57), ]
g1.m <- femke.nonused[sample(rownames(femke.nonused[femke.nonused$gender=="male" & femke.nonused$age.years<=40,]), 57), ]
g1 <- rbind(g1.f,g1.m)
g2.f <- femke.nonused[sample(rownames(femke.nonused[femke.nonused$gender=="female" & femke.nonused$age.years>40 & femke.nonused$age.years<=60,]), 48), ]
g2.m <- femke.nonused[sample(rownames(femke.nonused[femke.nonused$gender=="male" & femke.nonused$age.years>40 & femke.nonused$age.years<=60,]), 49), ]
g2 <- rbind(g2.f,g2.m)
g3.f <- femke.nonused[sample(rownames(femke.nonused[femke.nonused$gender=="female" & femke.nonused$age.years>60 & femke.nonused$age.years<80,]), 21), ]
g3.m <- femke.nonused[sample(rownames(femke.nonused[femke.nonused$gender=="male" & femke.nonused$age.years>60 & femke.nonused$age.years<80,]), 22), ]
g3 <- rbind(g3.f,g3.m)
#g4.f <- femke.nonused[sample(rownames(femke.nonused[femke.nonused$gender=="female" & femke.nonused$age.years>=80,]), 40), ]
#g4.m <- femke.nonused[sample(rownames(femke.nonused[femke.nonused$gender=="male" & femke.nonused$age.years>=80,]), 40), ]
#g4 <- rbind(g4.f,g4.m)
g <- rbind(g1,g2,g3)
gc <- g[,c("new_id","gender","age.years")]
colnames(gc) <- colnames(data.2015)
gt <- rbind(data.2015,gc)
hist(data.2015$Age[data.2015$Gender=="female"],main="Females histogram 2015",xlab="age",ylim=c(0,35))
hist(data.2015$Age[data.2015$Gender=="male"],main="Males histogram 2015",xlab="age",ylim=c(0,35))
hist(gt$Age[gt$Gender=="female"],main="Females histogram + requested 2016",xlab="age")
hist(gt$Age[gt$Gender=="male"],main="Males histogram + requested2016",xlab="age")
write.table(gc,file=paste(file.dir,"16Srequest2016.txt",sep="/"),col.names = T,sep="\t",quote = F)
write.table(gt,file=paste(file.dir,"16S2015plusRequested2016.txt",sep="/"),col.names = T,sep="\t",quote = F)

metag_miss2015 <-  subset(data.2015, !data.2015$ID%in%last49$Verschickungscode)
write.table(metag_miss2015,file=paste(file.dir,"metagenome_miss.txt",sep="/"),col.names = T,sep="\t",quote = F)


samples.unique <- unique(samples)
samples.duplicate <- as.data.frame(samples[duplicated(samples),])
write.table(samples.unique,'/home/torres/Documents/Projects/Metagenome/data/16s/samples.unique.txt',col.names = F,row.names = F,sep="\t",quote = F)
write.table(samples.duplicate,'/home/torres/Documents/Projects/Metagenome/data/16s/samples.duplicated.txt',col.names = F,row.names = F,sep="\t",quote = F)


