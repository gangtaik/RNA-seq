#parameter calling
args=commandArgs(trailingOnly = TRUE)
cat(args)
brca.dir=args[1]
#brca.dir="/data3/gangtl95/BRCA/"
geo=args[2]
#geo="GSE115591"
setwd(paste0(brca.dir,"breast.",geo,"/htseq/"))
getwd()
#setting path
#loc=zeta(MRC3)

#read file in directory
dir = c(getwd());dir
file_list = list.files(dir, pattern="htseq.txt");length(file_list)

for ( i in file_list){
  print(i)
  assign(i, read.table(i, sep="\t", stringsAsFactors = F, header = F))
}

# merge file
df = data.frame((get(file_list[1]))[1])
for ( i in file_list){
  t1 = get(i)
  colnames(t1) = c("V1", i)
  df = merge(df, t1, by="V1")
}
df = df[-c(1:5),]
# column name rename
colnames(df) = gsub(".htseq.txt", "", colnames(df))
rownames(df) = df$V1
df = df[,-1]
df[1:10,]

# write merging counting files
# b.f normalization
write.table(df, paste0(dir,"/breast.",geo,"_readcounts.txt"),sep = "\t", row.names = T, quote = F)
