library(ggplot2)
library(ramwas)
library(data.table)

get_res = function(rdir){
    data_lst = list()
    for (i in 1:22) {
        spath = paste0(rdir, i, '.sharepro.gxe.txt')
        data_lst[[i]] = fread(spath)
    }
    df = do.call(rbind, data_lst)
}

plot_man = function(nlogp, chr, pos, filename, ctf){
    pdf(file = filename, height = 4, width = 6)
    plt = manPlotPrepare(nlogp, chr, pos, ismlog10 = TRUE)
    manPlotFast(plt, colorSet = c("plum","darkblue"))
    abline(h = -log10(ctf),lty = 2,col = "red")
    dev.off()
}

plot_man_yl = function(nlogp, chr, pos, filename, ctf, yl){
    pdf(file = filename, height = 4, width = 6)
    plt = manPlotPrepare(nlogp, chr, pos, ismlog10 = TRUE)
    manPlotFast(plt, ylim = c(0, yl), colorSet = c("plum","darkblue"))
    abline(h = -log10(ctf),lty = 2,col = "red")
    dev.off()
}

get_clump = function(rdir){
    data_lst = list()
    for (i in 1:22) {
        spath = paste0(rdir, i, '.clumped')
        data_lst[[i]] = fread(spath)
    }
    df = do.call(rbind, data_lst)
}

get_cojo = function(rdir){
    data_lst = list()
    for (i in 1:22) {
        spath = paste0(rdir, i, '.jma.cojo')
        data_lst[[i]] = fread(spath)
    }
    df = do.call(rbind, data_lst)
}

WHRadjBMI = get_res('../dat/WHRadjBMI_sex_FM/WHRadjBMI_FM_')
WHRadjBMI_AL_clump = get_clump('../dat/WHRadjBMI_sex_FM/WHRadjBMI_AL_')
WHRadjBMI_AL_cojo = get_cojo('../dat/WHRadjBMI_sex_FM/WHRadjBMI_AL_')
WHRadjBMI$nlog10pGxE_0_1 = (pchisq((WHRadjBMI$BETA0-WHRadjBMI$BETA1)**2/(WHRadjBMI$SE0^2+WHRadjBMI$SE1^2), df=1, lower.tail = F,log.p = T)/-log(10))
WHRadjBMI$SharePro = 0
WHRadjBMI$Clump = 0
WHRadjBMI$COJO = 0
WHRadjBMI$SharePro[WHRadjBMI$cs!=0] = 1
WHRadjBMI$Clump[WHRadjBMI$SNP %in% WHRadjBMI_AL_clump$SNP] = 1
WHRadjBMI$COJO[WHRadjBMI$SNP %in% WHRadjBMI_AL_cojo$SNP] = 1

write.table(WHRadjBMI, file = '../doc/WHR_statistics.txt', quote=F, sep='\t', row.names=F, col.names=T)
plot_man(WHRadjBMI$nlog10pGxE_0_1, WHRadjBMI$CHR, WHRadjBMI$POS, '../fig/WHR_FM.pdf', 5e-8)
plot_man(WHRadjBMI$nlog10pGxE_0_1[WHRadjBMI$SharePro==1], WHRadjBMI$CHR[WHRadjBMI$SharePro==1], WHRadjBMI$POS[WHRadjBMI$SharePro==1], '../fig/WHR_FM_SharePro.pdf', 0.05/sum(WHRadjBMI$SharePro==1))
plot_man(WHRadjBMI$nlog10pGxE_0_1[WHRadjBMI$COJO==1], WHRadjBMI$CHR[WHRadjBMI$COJO==1], WHRadjBMI$POS[WHRadjBMI$COJO==1], '../fig/WHR_FM_COJO.pdf', 0.05/sum(WHRadjBMI$COJO==1))
plot_man(WHRadjBMI$nlog10pGxE_0_1[WHRadjBMI$Clump==1], WHRadjBMI$CHR[WHRadjBMI$Clump==1], WHRadjBMI$POS[WHRadjBMI$Clump==1], '../fig/WHR_FM_Clump.pdf', 0.05/sum(WHRadjBMI$Clump==1))


WHRadjBMI_INT = get_res('../dat/WHRadjBMI_sex_FM/WHRadjBMI_FM_INT_')
WHRadjBMI_AL_INT_clump = get_clump('../dat/WHRadjBMI_sex_FM/WHRadjBMI_AL_INT_')
WHRadjBMI_AL_INT_cojo = get_cojo('../dat/WHRadjBMI_sex_FM/WHRadjBMI_AL_INT_')
WHRadjBMI_INT$nlog10pGxE_0_1 = (pchisq((WHRadjBMI_INT$BETA0-WHRadjBMI_INT$BETA1)**2/(WHRadjBMI_INT$SE0^2+WHRadjBMI_INT$SE1^2), df=1, lower.tail = F,log.p = T)/-log(10))
WHRadjBMI_INT$SharePro = 0
WHRadjBMI_INT$Clump = 0
WHRadjBMI_INT$COJO = 0
WHRadjBMI_INT$SharePro[WHRadjBMI_INT$cs!=0] = 1
WHRadjBMI_INT$Clump[WHRadjBMI_INT$SNP %in% WHRadjBMI_AL_INT_clump$SNP] = 1
WHRadjBMI_INT$COJO[WHRadjBMI_INT$SNP %in% WHRadjBMI_AL_INT_cojo$SNP] = 1

write.table(WHRadjBMI_INT, file = '../doc/WHR_INT_statistics.txt', quote=F, sep='\t', row.names=F, col.names=T)
plot_man(WHRadjBMI_INT$nlog10pGxE_0_1, WHRadjBMI_INT$CHR, WHRadjBMI_INT$POS, '../fig/WHR_INT_FM.pdf', 5e-8)
plot_man(WHRadjBMI_INT$nlog10pGxE_0_1[WHRadjBMI_INT$SharePro==1], WHRadjBMI_INT$CHR[WHRadjBMI_INT$SharePro==1], WHRadjBMI_INT$POS[WHRadjBMI_INT$SharePro==1], '../fig/WHR_INT_FM_SharePro.pdf', 0.05/sum(WHRadjBMI_INT$SharePro==1))
plot_man(WHRadjBMI_INT$nlog10pGxE_0_1[WHRadjBMI_INT$COJO==1], WHRadjBMI_INT$CHR[WHRadjBMI_INT$COJO==1], WHRadjBMI_INT$POS[WHRadjBMI_INT$COJO==1], '../fig/WHR_INT_FM_COJO.pdf', 0.05/sum(WHRadjBMI_INT$COJO==1))
plot_man(WHRadjBMI_INT$nlog10pGxE_0_1[WHRadjBMI_INT$Clump==1], WHRadjBMI_INT$CHR[WHRadjBMI_INT$Clump==1], WHRadjBMI_INT$POS[WHRadjBMI_INT$Clump==1], '../fig/WHR_INT_FM_Clump.pdf', 0.05/sum(WHRadjBMI_INT$Clump==1))

FFR = get_res('/home/wmzh22/scratch/SharePro_gxe/dat/FFR_smoke_CEN/FFR_CEN_')
FFR$nlog10pGxE_0_1 = (pchisq((FFR$BETA0-FFR$BETA1)**2/(FFR$SE0^2+FFR$SE1^2), df=1, lower.tail = F,log.p = T)/-log(10))
FFR$nlog10pGxE_0_2 = (pchisq((FFR$BETA0-FFR$BETA2)**2/(FFR$SE0^2+FFR$SE2^2), df=1, lower.tail = F,log.p = T)/-log(10))
FFR$nlog10pGxE_1_2 = (pchisq((FFR$BETA1-FFR$BETA2)**2/(FFR$SE1^2+FFR$SE2^2), df=1, lower.tail = F,log.p = T)/-log(10))
FFR_AS_clump = get_clump('/home/wmzh22/scratch/SharePro_gxe/dat/FFR_smoke_CEN/FFR_AS_')
FFR_AS_cojo = get_cojo('/home/wmzh22/scratch/SharePro_gxe/dat/FFR_smoke_CEN/FFR_AS_')
FFR$SharePro = 0
FFR$Clump = 0
FFR$COJO = 0
FFR$SharePro[FFR$cs!=0] = 1
FFR$Clump[FFR$SNP %in% FFR_AS_clump$SNP] = 1
FFR$COJO[FFR$SNP %in% FFR_AS_cojo$SNP] = 1

write.table(FFR, file = '../doc/FFR_statistics.txt', quote=F, sep='\t', row.names=F, col.names=T)
plot_man_yl(FFR$nlog10pGxE_0_1, FFR$CHR, FFR$POS, '../fig/FFR_CE.pdf', 5e-8, 20)
plot_man_yl(FFR$nlog10pGxE_0_2, FFR$CHR, FFR$POS, '../fig/FFR_CN.pdf', 5e-8, 20)
plot_man_yl(FFR$nlog10pGxE_1_2, FFR$CHR, FFR$POS, '../fig/FFR_EN.pdf', 5e-8, 20)
plot_man_yl(FFR$nlog10pGxE_0_1[FFR$SharePro==1], FFR$CHR[FFR$SharePro==1], FFR$POS[FFR$SharePro==1], '../fig/FFR_CE_SharePro.pdf', 0.05/sum(FFR$SharePro==1)/3, 20)
plot_man_yl(FFR$nlog10pGxE_0_2[FFR$SharePro==1], FFR$CHR[FFR$SharePro==1], FFR$POS[FFR$SharePro==1], '../fig/FFR_CN_SharePro.pdf', 0.05/sum(FFR$SharePro==1)/3, 20)
plot_man_yl(FFR$nlog10pGxE_1_2[FFR$SharePro==1], FFR$CHR[FFR$SharePro==1], FFR$POS[FFR$SharePro==1], '../fig/FFR_EN_SharePro.pdf', 0.05/sum(FFR$SharePro==1)/3, 20)
plot_man_yl(FFR$nlog10pGxE_0_1[FFR$Clump==1], FFR$CHR[FFR$Clump==1], FFR$POS[FFR$Clump==1], '../fig/FFR_CE_Clump.pdf', 0.05/sum(FFR$Clump==1)/3, 20)
plot_man_yl(FFR$nlog10pGxE_0_2[FFR$Clump==1], FFR$CHR[FFR$Clump==1], FFR$POS[FFR$Clump==1], '../fig/FFR_CN_Clump.pdf', 0.05/sum(FFR$Clump==1)/3, 20)
plot_man_yl(FFR$nlog10pGxE_1_2[FFR$Clump==1], FFR$CHR[FFR$Clump==1], FFR$POS[FFR$Clump==1], '../fig/FFR_EN_Clump.pdf', 0.05/sum(FFR$Clump==1)/3, 20)
plot_man_yl(FFR$nlog10pGxE_0_1[FFR$COJO==1], FFR$CHR[FFR$COJO==1], FFR$POS[FFR$COJO==1], '../fig/FFR_CE_COJO.pdf', 0.05/sum(FFR$COJO==1)/3, 20)
plot_man_yl(FFR$nlog10pGxE_0_2[FFR$COJO==1], FFR$CHR[FFR$COJO==1], FFR$POS[FFR$COJO==1], '../fig/FFR_CN_COJO.pdf', 0.05/sum(FFR$COJO==1)/3, 20)
plot_man_yl(FFR$nlog10pGxE_1_2[FFR$COJO==1], FFR$CHR[FFR$COJO==1], FFR$POS[FFR$COJO==1], '../fig/FFR_EN_COJO.pdf', 0.05/sum(FFR$COJO==1)/3, 20)

FFR_INT = get_res('/home/wmzh22/scratch/SharePro_gxe/dat/FFR_smoke_CEN/FFR_CEN_INT_')
FFR_INT$nlog10pGxE_0_1 = (pchisq((FFR_INT$BETA0-FFR_INT$BETA1)**2/(FFR_INT$SE0^2+FFR_INT$SE1^2), df=1, lower.tail = F, log.p = T)/-log(10))
FFR_INT$nlog10pGxE_0_2 = (pchisq((FFR_INT$BETA0-FFR_INT$BETA2)**2/(FFR_INT$SE0^2+FFR_INT$SE2^2), df=1, lower.tail = F, log.p = T)/-log(10))
FFR_INT$nlog10pGxE_1_2 = (pchisq((FFR_INT$BETA1-FFR_INT$BETA2)**2/(FFR_INT$SE1^2+FFR_INT$SE2^2), df=1, lower.tail = F, log.p = T)/-log(10))
FFR_INT_AS_clump = get_clump('/home/wmzh22/scratch/SharePro_gxe/dat/FFR_smoke_CEN/FFR_AS_INT_')
FFR_INT_AS_cojo = get_cojo('/home/wmzh22/scratch/SharePro_gxe/dat/FFR_smoke_CEN/FFR_AS_INT_')
FFR_INT$SharePro = 0
FFR_INT$Clump = 0
FFR_INT$COJO = 0
FFR_INT$SharePro[FFR_INT$cs!=0] = 1
FFR_INT$Clump[FFR_INT$SNP %in% FFR_INT_AS_clump$SNP] = 1
FFR_INT$COJO[FFR_INT$SNP %in% FFR_INT_AS_cojo$SNP] = 1

write.table(FFR_INT, file = '../doc/FFR_INT_statistics.txt', quote=F, sep='\t', row.names=F, col.names=T)

plot_man_yl(FFR_INT$nlog10pGxE_0_1, FFR_INT$CHR, FFR_INT$POS, '../fig/FFR_INT_CE.pdf', 5e-8, 20)
plot_man_yl(FFR_INT$nlog10pGxE_0_2, FFR_INT$CHR, FFR_INT$POS, '../fig/FFR_INT_CN.pdf', 5e-8, 20)
plot_man_yl(FFR_INT$nlog10pGxE_1_2, FFR_INT$CHR, FFR_INT$POS, '../fig/FFR_INT_EN.pdf', 5e-8, 20)
plot_man_yl(FFR_INT$nlog10pGxE_0_1[FFR_INT$SharePro==1], FFR_INT$CHR[FFR_INT$SharePro==1], FFR_INT$POS[FFR_INT$SharePro==1], '../fig/FFR_INT_CE_SharePro.pdf', 0.05/sum(FFR_INT$SharePro==1)/3, 20)
plot_man_yl(FFR_INT$nlog10pGxE_0_2[FFR_INT$SharePro==1], FFR_INT$CHR[FFR_INT$SharePro==1], FFR_INT$POS[FFR_INT$SharePro==1], '../fig/FFR_INT_CN_SharePro.pdf', 0.05/sum(FFR_INT$SharePro==1)/3, 20)
plot_man_yl(FFR_INT$nlog10pGxE_1_2[FFR_INT$SharePro==1], FFR_INT$CHR[FFR_INT$SharePro==1], FFR_INT$POS[FFR_INT$SharePro==1], '../fig/FFR_INT_EN_SharePro.pdf', 0.05/sum(FFR_INT$SharePro==1)/3, 20)
plot_man_yl(FFR_INT$nlog10pGxE_0_1[FFR_INT$Clump==1], FFR_INT$CHR[FFR_INT$Clump==1], FFR_INT$POS[FFR_INT$Clump==1], '../fig/FFR_INT_CE_Clump.pdf', 0.05/sum(FFR_INT$Clump==1)/3, 20)
plot_man_yl(FFR_INT$nlog10pGxE_0_2[FFR_INT$Clump==1], FFR_INT$CHR[FFR_INT$Clump==1], FFR_INT$POS[FFR_INT$Clump==1], '../fig/FFR_INT_CN_Clump.pdf', 0.05/sum(FFR_INT$Clump==1)/3, 20)
plot_man_yl(FFR_INT$nlog10pGxE_1_2[FFR_INT$Clump==1], FFR_INT$CHR[FFR_INT$Clump==1], FFR_INT$POS[FFR_INT$Clump==1], '../fig/FFR_INT_EN_Clump.pdf', 0.05/sum(FFR_INT$Clump==1)/3, 20)
plot_man_yl(FFR_INT$nlog10pGxE_0_1[FFR_INT$COJO==1], FFR_INT$CHR[FFR_INT$COJO==1], FFR_INT$POS[FFR_INT$COJO==1], '../fig/FFR_INT_CE_COJO.pdf', 0.05/sum(FFR_INT$COJO==1)/3, 20)
plot_man_yl(FFR_INT$nlog10pGxE_0_2[FFR_INT$COJO==1], FFR_INT$CHR[FFR_INT$COJO==1], FFR_INT$POS[FFR_INT$COJO==1], '../fig/FFR_INT_CN_COJO.pdf', 0.05/sum(FFR_INT$COJO==1)/3, 20)
plot_man_yl(FFR_INT$nlog10pGxE_1_2[FFR_INT$COJO==1], FFR_INT$CHR[FFR_INT$COJO==1], FFR_INT$POS[FFR_INT$COJO==1], '../fig/FFR_INT_EN_COJO.pdf', 0.05/sum(FFR_INT$COJO==1)/3, 20)


FFR_sig = FFR[FFR$cs!=0 & (FFR$nlog10pdiff_0_1 > -log10(0.05/sum(FFR$cs!=0)/3) | FFR$nlog10pdiff_0_2 > -log10(0.05/sum(FFR$cs!=0)/3) | FFR$nlog10pdiff_1_2 > -log10(0.05/sum(FFR$cs!=0)/3)),]
gtf = as.data.frame(fread("../dat/gencode.v19.protein_coding.gtf"))
gtf$gene = gsub("\"", "", gsub(" gene_name ", "", unlist(lapply(strsplit(gtf$V9, ";"), function(x) x[5]))))

get_nearest_gene = function(gtf, chr, pos){
    subgtf = gtf[gtf$V1 == paste0("chr", chr),]
    intersectgene = subgtf[subgtf$V4 <= pos & subgtf$V5 >= pos,]
    if (nrow(intersectgene)!=0) {
        ngene = paste(intersectgene$gene, collapse=";")
    } else {
        subgtf$leftdist = abs(subgtf$V4 - pos)
        subgtf$rightdist = abs(subgtf$V5 - pos)
        subgtf$mindist = apply(subgtf[, c("leftdist","rightdist")], 1, function(x) min(x))
        if (min(subgtf$mindist) > 500000) {
            ngene = 'NA'
        } else {
            ngene = paste(subgtf[subgtf$mindist == min(subgtf$mindist),]$gene, collapse=";")
        }
    }
    return(ngene)
}
FFR_sig$NGene = mapply(get_nearest_gene, chr = FFR_sig$CHR, pos = FFR_sig$POS, MoreArgs = list(gtf = gtf))
pltdat = data.frame(Beta = c(FFR_sig$BETA0, FFR_sig$BETA1, FFR_sig$BETA2), SE = c(FFR_sig$SE0, FFR_sig$SE1, FFR_sig$SE2), population = rep(c("Current smokers","Previous smokers","Never smokers"), each = nrow(FFR_sig)))
pltdat$NGene = rep(FFR_sig$NGene, 3) 
pltdat$population = factor(pltdat$population, levels = c("Current smokers","Previous smokers","Never smokers"))
pltdat$Annotation = paste0(pltdat$NGene," (",paste0(FFR_sig$CHR,":",FFR_sig$POS,"_",FFR_sig$A1,"_",FFR_sig$A2),")")
pltdat$Annotation = factor(pltdat$Annotation, levels = rev(unique(pltdat$Annotation)))

ggplot(pltdat, aes(x = Beta, y = Annotation, color = population)) +
  geom_vline(xintercept = 0, lty = 2, color = "darkgrey") +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_errorbarh(aes(xmin = Beta - 1.96 * SE,
                     xmax = Beta + 1.96 * SE),height = 0,position = position_dodge(width = 0.3)) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.position = "right") +
  scale_color_manual(values = c("darkred", "orange", "royalblue")) +
  labs(color = "") +
  xlab(expression(hat(beta)[FFR])) +
  ylab("Nearest Gene") -> plt
ggsave("../fig/FFR_smoking_gene.pdf", plt, height = 6, width = 8)