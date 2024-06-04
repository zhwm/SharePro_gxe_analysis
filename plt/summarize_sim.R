library(PRROC)
library(ggplot2)
library(data.table)
library(ramwas)

wdir = "../sim/"
odir = "../doc/"
fdir = "../fig/"
setwd(wdir)

folderlist = c('3_0.05_-0.05/', '3_0.05_-0.02/', '3_0.05_0/', '3_0.05_0.02/', '3_0.05_0.05/', '2_0.05_-0.05/', '2_0.05_-0.02/', '2_0.05_0/', '2_0.05_0.02/', '2_0.05_0.05/', '1_0.05_-0.05/', '1_0.05_-0.02/', '1_0.05_0/', '1_0.05_0.02/', '1_0.05_0.05/')
Locilist = c("Locus1/", "Locus2/", "Locus3/")

get_auprc_vec = function(pip,gt){
    cv = pr.curve(pip[gt],pip[!gt])
    return(cv$auc.integral)
}

get_auroc_vec = function(pip,gt){
    cv = roc.curve(pip[gt],pip[!gt])
    return(cv$auc)
}

auprc_pip_vec = c()
auroc_pip_vec = c()
for (folder in folderlist) {
    params = unlist(strsplit(gsub('/','',folder),'_'))
    for (rt in c('CL','CM','CS')) {
        for (lo in Locilist) {
            SApip = c()
            Supip = c()
            SHpip = c()
            pval1 = c()
            pval2 = c()
            GTv = c()
            csnp = read.table(paste0(wdir,lo,folder,'csnp.txt'),header=F)
            cat(paste0(lo,folder,rt),collapse="\n")
            for (x in 1:50) {
                gem = read.table(paste0(wdir,lo,folder,rt,x), header=T)
                rownames(gem) = gem$SNPID
                SH = fread(paste0(wdir,lo,folder,rt,x,'.sharepro.gxe.txt'))
                SA = fread(paste0(wdir,lo,folder,rt,x,'.sharepro.combine.txt'))
                SU = fread(paste0(wdir,lo,folder,rt,x,'.susie.txt'))
                SApip = c(SApip, SA$PIP)
                Supip = c(Supip, SU$PIP)
                SHpip = c(SHpip, SH$PIP)
                pval1 = c(pval1, -gem$P_Value_Marginal)
                pval2 = c(pval2, -gem$P_Value_Joint)
                GTv = c(GTv,gem$SNPID %in% csnp[x,])
            }
            pips = list(SApip,Supip,SHpip,pval1,pval2)
            auprc_pip_vec = c(auprc_pip_vec, unlist(lapply(pips, function(y){get_auprc_vec(y,GTv)})))
            auprc_pip_vec = c(auprc_pip_vec, mean(GTv))
            auroc_pip_vec = c(auroc_pip_vec, unlist(lapply(pips, function(y){get_auroc_vec(y,GTv)})))
            auroc_pip_vec = c(auroc_pip_vec, 0.5)
        }
    }
}

auprc = data.frame(AUPRC = auprc_pip_vec, Method = c("SparsePro", "SuSiE", "SharePro", "GEM-1df", "GEM-2df", "Baseline"), Loci=rep(Locilist, each=6), Ratio = rep(c('CL','CM','CS'), each=6*length(Locilist)), Folder = rep(folderlist, each = 6*length(Locilist)*3))
auprc$K = unlist(lapply(strsplit(auprc$Folder,'_'),function(x){x[1]}))
auprc$Ne = 25000
auprc$Nu = factor(unlist(lapply(auprc$Ratio, function(x){strsplit(x,2,2)})), labels = c('25000','10000','5000'))
auprc$Be = unlist(lapply(strsplit(auprc$Folder,'_'),function(x){x[2]}))
auprc$Bu = unlist(lapply(strsplit(auprc$Folder,'_'),function(x){gsub('/','',x[3])}))
auprc$Method = factor(auprc$Method, levels = c("GEM-1df","GEM-2df","SuSiE","SparsePro","SharePro", "Baseline"))
write.table(auprc, file = '../doc/sim_auprc.txt', quote=F, sep='\t', row.names=F, col.names=T)

auroc = data.frame(AUROC = auroc_pip_vec, Method = c("SparsePro", "SuSiE", "SharePro", "GEM-1df", "GEM-2df", "Baseline"), Loci=rep(Locilist, each=6), Ratio = rep(c('CL','CM','CS'), each=6*length(Locilist)), Folder = rep(folderlist, each = 6*length(Locilist)*3))
auroc$K = unlist(lapply(strsplit(auroc$Folder,'_'),function(x){x[1]}))
auroc$Ne = 25000
auroc$Nu = factor(unlist(lapply(auroc$Ratio, function(x){strsplit(x,2,2)})), labels = c('25000','10000','5000'))
auroc$Be = unlist(lapply(strsplit(auroc$Folder,'_'),function(x){x[2]}))
auroc$Bu = unlist(lapply(strsplit(auroc$Folder,'_'),function(x){gsub('/','',x[3])}))
auroc$Method = factor(auroc$Method, levels = c("GEM-1df","GEM-2df","SuSiE","SparsePro","SharePro", "Baseline"))
write.table(auroc, file = '../doc/sim_auroc.txt', quote=F, sep='\t', row.names=F, col.names=T)
