library(susieR)
library(data.table)
args = commandArgs(trailingOnly=T)
file = args[1]
ld = args[2]
ss = read.table(file, header=T)
ld = fread(ld)
model = susie_rss(ss$Beta_Marginal/ss$SE_Beta_Marginal, as.matrix(ld))
res = data.frame(SNP = ss$SNPID)
res$PIP = model$pip
write.table(res,paste0(file,".susie.txt"), col.names=T,row.names=F,sep="\t",quote=F)
