library(data.table)

wdir = "../sim/"
odir = "../doc/"
fdir = "../fig/"
setwd(wdir)

pfolderlist = c('3_0.05_-0.05/', '3_0.05_-0.02/', '3_0.05_0/', '3_0.05_0.02/', '2_0.05_-0.05/', '2_0.05_-0.02/', '2_0.05_0/', '2_0.05_0.02/', '1_0.05_-0.05/', '1_0.05_-0.02/', '1_0.05_0/', '1_0.05_0.02/')
Locilist = c("Locus1/", "Locus2/", "Locus3/")

get_power = function(pval, Sv, GTv, ctf, adj = "BH"){
    pvalselected = pval[Sv==1]
    pvalselected = p.adjust(pvalselected, method = adj)
    totalp = sum(GTv)
    truep = sum(pvalselected < ctf & GTv[Sv==1])
    return(truep / totalp)
}

get_fdr = function(pval, Sv, GTv, ctf, adj = "BH"){
    pvalselected = pval[Sv==1]
    pvalselected = p.adjust(pvalselected, method = adj)
    falsep = sum(pvalselected < ctf & !GTv[Sv==1])
    predp = sum(pvalselected < ctf)
    if (predp > 0) {
        fdr = falsep / predp
    } else {
        fdr = 0.0
    }
    return(fdr)
}

ctf = 0.01
fdr_pval_vec = c()
pwr_pval_vec = c()
for (folder in pfolderlist) {
    for (rt in c('CL','CM','CS')) {
        tp = 0
        predp = 0
        allp = 0
        gxegem = c()
        gxeadj = c()
        gxecojo = c()
        gxeclump = c()
        gxesharepro = c()
        GTv = c()
        for (lo in Locilist) {
        cat(paste0(lo,folder,rt),collapse="\n")
            csnp = read.table(paste0(wdir,lo,folder,'csnp.txt'),header=F)
            for (x in 1:50) {
                gem = read.table(paste0(wdir,lo,folder,rt,x), header=T)
                rownames(gem) = gem$SNPID
                gem$cojop = 0
                gem$clumpp = 0
                gem$shp = 0
                gem$adj = 0
                gem$adj[p.adjust(gem$P_Value_Interaction, 'fdr')<ctf] = 1
                if (paste0(rt,x,'.jma.cojo') %in% list.files(paste0(wdir,lo,folder))) {
                    cojo = fread(paste0(wdir,lo,folder,rt,x,'.jma.cojo'))
                    gem[cojo$SNP,'cojop'] = 1
                }
                if (paste0(rt,x,'.clumped') %in% list.files(paste0(wdir,lo,folder))) {
                    clump = fread(paste0(wdir,lo,folder,rt,x,'.clumped'))
                    gem[clump$SNP,'clumpp'] = 1
                }
                SH = fread(paste0(wdir,lo,folder,rt,x,'.sharepro.gxe.txt'))
                gem[SH$cs!=0,'shp'] = 1
                gxegem = c(gxegem, gem$P_Value_Interaction)
                gxeadj = c(gxeadj, gem$adj)
                gxecojo = c(gxecojo, gem$cojop)
                gxeclump = c(gxeclump, gem$clumpp)
                gxesharepro = c(gxesharepro, gem$shp)
                GTv = c(GTv,gem$SNPID %in% csnp[x,])
                allp = allp + length(csnp[x, ])
                cond = SH$cs!=0 & gem$P_Value_Interaction<ctf
                if (sum(cond)>0) {
                    for (e in strsplit(SH$cs_variants[cond], '/')) {
                        if (length(intersect(e, csnp[x,]))>0) {
                            tp = tp + 1
                        }
                    }
                }
                predp = predp + sum(cond)
            }
        }
        idxs = list(gxeadj,gxecojo,gxeclump,gxesharepro)
        fdr_pval_vec = c(fdr_pval_vec, 1-tp/predp)
        fdr_pval_vec = c(fdr_pval_vec, unlist(lapply(idxs, function(x){get_fdr(gxegem, x, GTv, ctf)})))
        pwr_pval_vec = c(pwr_pval_vec, tp/allp)
        pwr_pval_vec = c(pwr_pval_vec, unlist(lapply(idxs, function(x){get_power(gxegem, x, GTv, ctf)})))
    }
}

df_pwr_fdr = data.frame(Power = pwr_pval_vec, FDR = fdr_pval_vec, Method = c("SharePro(effect)", "GEM", "COJO", "Clump", "SharePro(variant)"), Ratio = rep(c('CL','CM','CS'), each=5), Folder = rep(pfolderlist, each = 5*3))
df_pwr_fdr$K = unlist(lapply(strsplit(df_pwr_fdr$Folder,'_'),function(x){x[1]}))
df_pwr_fdr$Ne = 25000
df_pwr_fdr$Nu = factor(unlist(lapply(df_pwr_fdr$Ratio, function(x){strsplit(x,2,2)})), labels = c('25000','10000','5000'))
df_pwr_fdr$Be = unlist(lapply(strsplit(df_pwr_fdr$Folder,'_'),function(x){x[2]}))
df_pwr_fdr$Bu = unlist(lapply(strsplit(df_pwr_fdr$Folder,'_'),function(x){gsub('/','',x[3])}))
write.table(df_pwr_fdr, file = '../doc/sim_pwr_fdr.txt', quote=F, sep='\t', row.names=F, col.names=T)
