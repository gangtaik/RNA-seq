#parameter calling
args=commandArgs(trailingOnly = TRUE)
cat(args)
brca.dir=args[1]
#brca.dir="/data/Epsilon_data12/gangtl95/BRCA/"
geo=args[2]
#geo="GSE150997"
setwd(paste0(brca.dir,"breast.",geo,"/hisat2.GRCh38/strandness"))
#rm(list=ls())
#Flag Check
.libPaths()
if(!require(dplyr)){install.packages("dplyr")};library(dplyr)
if(!require(ggplot2)){BiocManager::install("ggplot2")};library(ggplot2)
if(!require(gridExtra)){install.packages("gridExtra")};library(gridExtra)
options(scipen=100)

actb.list=dir(path =getwd(),pattern = "sorted.ACTB.txt")
gapdf.list=dir(path =getwd(),pattern = "sorted.GAPDH.txt")
numa1.list=dir(path =getwd(),pattern = "sorted.NUMA1.txt")
mt.list=dir(path =getwd(),pattern = "_pos.txt")


df.ACTB=c()
for (i in 1:length(actb.list)){
  a=sub(".txt",replacement = "",actb.list[i])
  print (a)
  assign(a,as.data.frame(read.table(actb.list[i], sep="\t",stringsAsFactors = F,header = F)))
  df.smp=data.frame(rep(a,nrow(get(a))),get(a))
  colnames(df.smp)=c("Sample","Num","Flag.ACTB")
  df.ACTB=rbind(df.ACTB,df.smp)
  
}
df.ACTB$Flag.ACTB = as.factor(df.ACTB$Flag.ACTB)

df.GAPDH=c()
for (i in 1:length(gapdf.list)){
  a=sub(".txt",replacement = "",gapdf.list[i])
  print (a)
  assign(a,as.data.frame(read.table(gapdf.list[i], sep="\t",stringsAsFactors = F,header = F)))
  df.smp=data.frame(rep(a,nrow(get(a))),get(a))
  colnames(df.smp)=c("Sample","Num","Flag.GAPDH")
  df.GAPDH=rbind(df.GAPDH,df.smp)
  
}
df.GAPDH$Flag.GAPDH = as.factor(df.GAPDH$Flag.GAPDH)

df.NUMA1=c()
for (i in 1:length(numa1.list)){
  a=sub(".txt",replacement = "",numa1.list[i])
  print (a)
  assign(a,as.data.frame(read.table(numa1.list[i], sep="\t",stringsAsFactors = F,header = F)))
  df.smp=data.frame(rep(a,nrow(get(a))),get(a))
  colnames(df.smp)=c("Sample","Num","Flag.NUMA1")
  df.NUMA1=rbind(df.NUMA1,df.smp)
  
}
df.NUMA1$Flag.NUMA1 = as.factor(df.NUMA1$Flag.NUMA1)

y.max=ceiling(log(max(df.ACTB$Num,df.NUMA1$Num,df.GAPDH$Num)+1,10))

png(filename=paste0(brca.dir,"breast.",geo,"/results/",geo,"_flag.png"),width = 500, height = 500, type = "cairo")

p= ggplot(df.ACTB, aes(x=Flag.ACTB, y=log(Num+1,base = 10)))+
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2),cex=0.6) +
  theme_bw() +
  scale_y_continuous(limits=c(0,y.max),breaks=c(0:y.max))

q = ggplot(df.GAPDH, aes(x=Flag.GAPDH, y=log(Num+1,base = 10))) + 
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2),cex=0.6) +
  theme_bw()+
  scale_y_continuous(limits=c(0,y.max),breaks=c(0:y.max))


r = ggplot(df.NUMA1, aes(x=Flag.NUMA1, y=log(Num+1,base = 10))) + 
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2),cex=0.6) +
  theme_bw()+
  scale_y_continuous(limits=c(0,y.max),breaks=c(0:y.max))

print(grid.arrange(p, q, r, ncol=2, nrow = 2))

invisible(dev.off())





#alignment 1 time reads 갯수 구하기 (라인 갯수)와 chr 비율 구하기
chr.list=c()
for (i in 1:length(mt.list)){
  a=sub("_pos.txt",replacement = "",mt.list[i])
  assign(a,as.data.frame(read.table(mt.list[i], sep="\t",stringsAsFactors = F,header = F)))
  df.smp=data.frame(paste("chr",get(a)[,2],sep=""),get(a)[,1]/sum(get(a)[,1]))
  get(a)[1:26,]
  chr.list=append(chr.list,get(a)[,2])
  chr.list=unique(chr.list)
  
}
chr.list=as.data.frame(chr.list)
colnames(chr.list)="Pos"
df.frq=chr.list

df.MT.percent=c()
for (i in 1:length(mt.list)){
  #i=1
  a=sub("_pos.txt",replacement = "",mt.list[i])
  assign(a,as.data.frame(read.table(mt.list[i], sep="\t",stringsAsFactors = F,header = F)))
  
  df.per=c()
  for (j in 1:nrow(get(a))){
    df.per=append(df.per,round(get(a)[j,1]/sum(get(a)[1:26,1])*100,4))
  }
  df.al=cbind(get(a),df.per)
  df.al=df.al[,c(2,3)]
  colnames(df.al)=c("Pos",a)
  df.frq=left_join(df.frq,df.al,by="Pos")
  
  
  #MT_percent
  df.smp=data.frame(paste("chr",get(a)[,2],sep=""),get(a)[,1]/sum(get(a)[,1]))
  colnames(df.smp)=c("Pos","Percent")
  mt=as.data.frame(cbind(sub("_pos.txt",replacement = "",mt.list[i]), df.smp$Percent[df.smp$Pos=="chrMT"]))
  colnames(mt)=c("Sample","MT.percent")
  df.MT.percent=rbind(df.MT.percent,mt)
  colnames(df.MT.percent)=colnames(mt)
}
df.frq[is.na(df.frq)]=0

write.table(df.MT.percent, file = paste0(brca.dir,"breast.",geo,"/results/",geo,"_MT_per.txt"),sep = "\t", row.names = F, quote = F)
write.table(df.frq, file = paste0(brca.dir,"breast.",geo,"/results/",geo,"_chr_per.txt"),sep = "\t", row.names = F, quote = F)

