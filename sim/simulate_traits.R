args = commandArgs(trailingOnly=TRUE) #LOCI KC KS hc hs
set.seed(2023)
library(data.table)

LO = args[1]
KC = as.numeric(args[2])
KS = as.numeric(args[3])
hf = as.numeric(args[4])
hl = as.numeric(args[5])
hs = as.numeric(args[6])
betaf = betal = rep(hs,(KC+2*KS))
betaf[1:(KC+KS)] = hf
betal[(KS+1):(KC+KS*2)] = hl
betadiff = betal - betaf
ite=50
print(betaf)
print(betal)

f = fread(paste0(LO,".raw"))
cond = rep(c(0,1),each=dim(f)[1]/2)

allcsnp = list()
ally = list(f$FID,f$IID,cond)

for (i in 1:ite){
  csnp = colnames(f)[sample(7:ncol(f),(KC+2*KS))]
  geno = scale(f[,csnp,with=F])
  genocond = geno*cond
  yhat = as.matrix(geno) %*% betaf + as.matrix(genocond) %*% betadiff
  y = yhat + rnorm(length(yhat),0,sqrt(1-var(yhat,na.rm = T)))
  ally[[i+3]]=as.vector(y)
  ally[[i+3]][is.na(ally[[i+3]])] <- 'NA'
  allcsnp[[i]]=unlist(lapply(strsplit(csnp,'_'),function(x){x[1]}))
}

setDT(ally)
setnames(ally, c(c('FID','IID','cond'),paste0('trait',seq(1,ite))))
fwrite(ally[1:30000],file='S.phen',sep='\t',quote=F,col.names=T) #Pheno file
fwrite(ally[1:35000],file='M.phen',sep='\t',quote=F,col.names=T) #Pheno file
fwrite(ally,file='L.phen',sep='\t',quote=F,col.names=T) #Pheno file
write.table(do.call(rbind,allcsnp),file='csnp.txt',row.names = F,col.names = F,sep='\t',quote=F) #True causal variants
