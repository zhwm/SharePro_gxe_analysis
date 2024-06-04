library(PRROC)
library(ggplot2)
library(data.table)
library(ramwas)

wdir = "../sim/"
odir = "../doc/"
fdir = "../fig/"
setwd(wdir)

lo = 'Locus1/'
folder = '1_0.05_0.02/'
rt = 'CL'
x = 21
bim = read.table(paste0(wdir,lo,gsub('/','.bim',lo)))
ld = read.table(paste0(wdir,lo,gsub('/','.ld',lo)))
csnp = read.table(paste0(wdir,lo,folder,'csnp.txt'),header=F)
gem = fread(paste0(wdir,lo,folder,rt,x))
plink = fread(paste0(wdir,lo,folder,rt,x,'.qassoc.gxe'))
plink$pval1 = pchisq((plink$BETA1/plink$SE1)^2,1,lower.tail = F)
plink$pval2 = pchisq((plink$BETA2/plink$SE2)^2,1,lower.tail = F)
SH = fread(paste0(wdir,lo,folder,rt,x,'.sharepro.gxe.txt'))
SA = fread(paste0(wdir,lo,folder,rt,x,'.sharepro.combine.txt'))

cslist = strsplit(SH$cs_variants[SH$cs!=0],'/')
gem$shp = NA
for (e in cslist) {
    gem$shp[gem$SNPID %in% unlist(e)] = gem$P_Value_Interaction[gem$SNPID==e[1]]
}

pvalplt = data.frame(chr = rep(gem$CHR,5), pos = rep(gem$POS,5), 
                    P = c(plink$pval1, 
                    plink$pval2,
                    gem$P_Value_Interaction,
                    p.adjust(gem$P_Value_Interaction,'fdr'),
                    gem$shp),
                    LD = rep(ld[,which(gem$SNPID==csnp[x,])]^2,5),
                    Adjustment = rep(c("Exposed association test",
                                        "Unexposed association test",
                                        "GEM interaction test",
                                        "GEM interaction test (FDR-adjusted)",
                                        "SharePro interaction test"),each=nrow(gem)))

ggplot(pvalplt[pvalplt$Adjustment=="Exposed association test" | pvalplt$Adjustment=="Unexposed association test",], aes(x = pos / 1000000, y = -log10(P), color=LD)) +
    facet_wrap( .~Adjustment,nrow = 5,strip.position = "top") +
    geom_point() +
    theme_classic() +
    theme(legend.position = "none",
        axis.text = element_text(size=13),
        axis.title = element_text(size=14),
        #      strip.background = element_blank(),
        strip.text = element_text(size = 13)) +
    xlab(paste0("Chr",pvalplt$chr[1]," (Mb)")) +
    ylab(expression(-log[10](p-value))) +
    ylim(0,20) +
    scale_color_stepsn(n.breaks = 6, colours = c("darkblue","blue","green","orange","red")) +
    labs(color = expression(r^2)) -> plt
ggsave(paste0(fdir, "example_E_UE_association.pdf"), plt, height = 5, width = 4)

ggplot(pvalplt[pvalplt$Adjustment=="GEM interaction test" | pvalplt$Adjustment=="GEM interaction test (FDR-adjusted)",], aes(x = pos / 1000000, y = -log10(P), color=LD)) +
    facet_wrap( .~Adjustment,nrow = 5,strip.position = "top") +
    geom_point() +
    theme_classic() +
    theme(legend.position = "none",
        axis.text = element_text(size=13),
        axis.title = element_text(size=14),
        #      strip.background = element_blank(),
        strip.text = element_text(size = 13)) +
    xlab(paste0("Chr",pvalplt$chr[1]," (Mb)")) +
    ylab(expression(-log[10](p-value))) +
    ylim(0,5) +
    scale_color_stepsn(n.breaks = 6, colours = c("darkblue","blue","green","orange","red")) +
    labs(color = expression(r^2)) -> plt
ggsave(paste0(fdir, "example_GEM_GxE.pdf"),plt,height = 5, width = 4)

ggplot(pvalplt[pvalplt$Adjustment=="SharePro interaction test",], aes(x = pos / 1000000, y = -log10(P), color=LD)) +
    facet_wrap( .~Adjustment,nrow = 5,strip.position = "top") +
    geom_point() +
    theme_classic() +
    theme(legend.position = "none",
        axis.text = element_text(size=13),
        axis.title = element_text(size=14),
        strip.text = element_text(size = 13)) +
    xlab(paste0("Chr",pvalplt$chr[1]," (Mb)")) +
    ylab(expression(-log[10](p-value))) +
    ylim(0,5) +
    scale_color_stepsn(n.breaks = 6, colours = c("darkblue","blue","green","orange","red")) +
    labs(color = expression(r^2)) -> plt
ggsave(paste0(fdir, "example_SharePro_GxE.pdf"),plt,height = 3, width = 4)

pipplt = data.frame(chr = rep(gem$CHR,1),
                    pos = rep(gem$POS,1),
                    P = c(SH$PIP),
                    LD = rep(ld[,which(gem$SNPID==csnp[x,])]^2,1),
                    Adjustment = rep(c("SharePro fine-mapping"),each=nrow(gem)))

ggplot(pipplt[pipplt$Adjustment=="SharePro fine-mapping",], aes(x = pos / 1000000, y = P, color=LD)) +
  facet_wrap( .~Adjustment,nrow = 5,strip.position = "top") +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size=13),
        axis.title = element_text(size=14),
        # strip.background = element_blank(),
        strip.text = element_text(size = 13)) +
  xlab(paste0("Chr",pipplt$chr[1]," (Mb)")) +
  ylab("PIP") + ylim(0,1) +
  scale_color_stepsn(n.breaks = 6, colours = c("darkblue","blue","green","orange","red")) +
  labs(color = expression(r^2)) -> plt
ggsave(paste0(fdir, "example_SharePro_finemap.pdf"),plt,height = 3, width = 4)

auprc = read.table('../doc/sim_auprc.txt', sep='\t', header=T)
auprc_formatted = data.frame(Number_of_causal_variants=auprc$K, 
                            Sample_size_in_exposed=auprc$Ne,
                            Sample_size_in_unexposed=auprc$Nu,
                            Effect_size_in_exposed=auprc$Be,
                            Effect_size_in_unexposed=auprc$Bu,
                            Loci=gsub('/', '', auprc$Loci),
                            Method=auprc$Method,
                            AUPRC=round(auprc$AUPRC,4))
write.table(auprc_formatted, file = '../doc/sim_auprc_formatted.txt', quote=F, sep='\t', row.names=F, col.names=T)

auprc$Method = factor(auprc$Method, levels = c("GEM-1df","GEM-2df","SuSiE","SparsePro","SharePro", "Baseline"))
auprc$Ne = paste0("N[e]:", auprc$Ne)
auprc$Nu = paste0("N[u]:", auprc$Nu)
auprc$Nu = factor(auprc$Nu, levels = c("N[u]:25000", "N[u]:10000", "N[u]:5000"))
auprc$Be = paste0("beta[e]:", auprc$Be)
auprc$Bu = paste0("beta[u]:", auprc$Bu)
auprc$Bu = factor(auprc$Bu, levels = c("beta[u]:-0.05", "beta[u]:-0.02", "beta[u]:0", "beta[u]:0.02", "beta[u]:0.05"))


ggplot(auprc[(auprc$K=='3') & (auprc$Method!='Baseline'),], aes(x=Method, y=AUPRC, color=Method)) + 
    facet_grid(Be+Bu~Ne+Nu,labeller = label_parsed) +
    geom_boxplot() + 
    geom_point() +
    theme_bw() +
    xlab("") + 
    ylab("AUPRC") +
    theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        title = element_text(size = 14),
        legend.position = "top") +
    scale_color_manual(values = c("orange","brown","lightblue", "royalblue","darkblue")) + 
    labs(color = "") -> plt
ggsave(paste0(fdir,"sim_auprc_3.pdf"),plt,height = 10, width = 10)

ggplot(auprc[(auprc$K=='2') & (auprc$Method!='Baseline'),], aes(x=Method, y=AUPRC, color=Method)) + 
    facet_grid(Be+Bu~Ne+Nu,labeller = label_parsed) +
    geom_boxplot() + 
    geom_point() +
    theme_bw() +
    xlab("") + 
    ylab("AUPRC") +
    theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        title = element_text(size = 14),
        legend.position = "top") +
    scale_color_manual(values = c("orange","brown","lightblue", "royalblue","darkblue")) + 
    labs(color = "") -> plt
ggsave(paste0(fdir,"sim_auprc_2.pdf"),plt,height = 10, width = 10)

ggplot(auprc[(auprc$K=='1') & (auprc$Method!='Baseline'),], aes(x=Method, y=AUPRC, color=Method)) + 
    facet_grid(Be+Bu~Ne+Nu,labeller = label_parsed) +
    geom_boxplot() + 
    geom_point() +
    theme_bw() +
    xlab("") + 
    ylab("AUPRC") +
    theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        title = element_text(size = 14),
        legend.position = "top") +
    scale_color_manual(values = c("orange","brown","lightblue", "royalblue","darkblue")) + 
    labs(color = "") -> plt
ggsave(paste0(fdir,"sim_auprc_1.pdf"),plt,height = 10, width = 10)

auroc = read.table('../doc/sim_auroc.txt', sep='\t', header=T)
auroc_formatted = data.frame(Number_of_causal_variants=auroc$K, 
                            Sample_size_in_exposed=auroc$Ne,
                            Sample_size_in_unexposed=auroc$Nu,
                            Effect_size_in_exposed=auroc$Be,
                            Effect_size_in_unexposed=auroc$Bu,
                            Loci=gsub('/', '', auroc$Loci),
                            Method=auroc$Method,
                            AUROC=round(auroc$AUROC,4))
write.table(auroc_formatted, file = '../doc/sim_auroc_formatted.txt', quote=F, sep='\t', row.names=F, col.names=T)

auroc$Method = factor(auroc$Method, levels = c("GEM-1df","GEM-2df","SuSiE","SparsePro","SharePro", "Baseline"))
auroc$Ne = paste0("N[e]:", auroc$Ne)
auroc$Nu = paste0("N[u]:", auroc$Nu)
auroc$Nu = factor(auroc$Nu, levels = c("N[u]:25000", "N[u]:10000", "N[u]:5000"))
auroc$Be = paste0("beta[e]:", auroc$Be)
auroc$Bu = paste0("beta[u]:", auroc$Bu)
auroc$Bu = factor(auroc$Bu, levels = c("beta[u]:-0.05", "beta[u]:-0.02", "beta[u]:0", "beta[u]:0.02", "beta[u]:0.05"))

ggplot(auroc[(auroc$K=='3') & (auroc$Method!='Baseline'),], aes(x=Method, y=AUROC, color=Method)) + 
    facet_grid(Be+Bu~Ne+Nu,labeller = label_parsed) +
    geom_boxplot() + 
    geom_point() +
    theme_bw() +
    xlab("") + 
    ylab("AUPRC") +
    theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        title = element_text(size = 14),
        legend.position = "top") +
    scale_color_manual(values = c("orange","brown","lightblue", "royalblue","darkblue")) + 
    labs(color = "") -> plt
ggsave(paste0(fdir,"sim_auroc_3.pdf"),plt,height = 10, width = 10)

ggplot(auroc[(auroc$K=='2') & (auroc$Method!='Baseline'),], aes(x=Method, y=AUROC, color=Method)) + 
    facet_grid(Be+Bu~Ne+Nu,labeller = label_parsed) +
    geom_boxplot() + 
    geom_point() +
    theme_bw() +
    xlab("") + 
    ylab("AUPRC") +
    theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        title = element_text(size = 14),
        legend.position = "top") +
    scale_color_manual(values = c("orange","brown","lightblue", "royalblue","darkblue")) + 
    labs(color = "") -> plt
ggsave(paste0(fdir,"sim_auroc_2.pdf"),plt,height = 10, width = 10)

ggplot(auroc[(auroc$K=='1') & (auroc$Method!='Baseline'),], aes(x=Method, y=AUROC, color=Method)) + 
    facet_grid(Be+Bu~Ne+Nu,labeller = label_parsed) +
    geom_boxplot() + 
    geom_point() +
    theme_bw() +
    xlab("") + 
    ylab("AUPRC") +
    theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        title = element_text(size = 14),
        legend.position = "top") +
    scale_color_manual(values = c("orange","brown","lightblue", "royalblue","darkblue")) + 
    labs(color = "") -> plt
ggsave(paste0(fdir,"sim_auroc_1.pdf"),plt,height = 10, width = 10)

df_pwr_fdr = read.table('../doc/sim_pwr_fdr.txt', sep='\t', header=T)

df_pwr_fdr_formatted = data.frame(Number_of_causal_variants=df_pwr_fdr$K, 
                            Sample_size_in_exposed=df_pwr_fdr$Ne,
                            Sample_size_in_unexposed=df_pwr_fdr$Nu,
                            Effect_size_in_exposed=df_pwr_fdr$Be,
                            Effect_size_in_unexposed=df_pwr_fdr$Bu,
                            Method=df_pwr_fdr$Method,
                            Power=round(df_pwr_fdr$Power,4),
                            FDR=round(df_pwr_fdr$FDR,4))

write.table(df_pwr_fdr_formatted, file = '../doc/sim_pwr_fdr_formatted.txt', quote=F, sep='\t', row.names=F, col.names=T)

df_pwr_fdr$Method = factor(df_pwr_fdr$Method, levels = c("GEM", "Clump", "COJO", "SharePro(variant)", "SharePro(effect)"))
df_pwr_fdr$Ne = paste0("N[e]:", df_pwr_fdr$Ne)
df_pwr_fdr$Nu = paste0("N[u]:", df_pwr_fdr$Nu)
df_pwr_fdr$Nu = factor(df_pwr_fdr$Nu, levels = c("N[u]:25000", "N[u]:10000", "N[u]:5000"))
df_pwr_fdr$Be = paste0("beta[e]:", df_pwr_fdr$Be)
df_pwr_fdr$Bu = paste0("beta[u]:", df_pwr_fdr$Bu)
df_pwr_fdr$Bu = factor(df_pwr_fdr$Bu, levels = c("beta[u]:-0.05", "beta[u]:-0.02", "beta[u]:0", "beta[u]:0.02", "beta[u]:0.05"))

ggplot(df_pwr_fdr[df_pwr_fdr$K==3,], aes(x=Power, y=FDR, color=Method, shape=Method)) + 
    geom_point(size=3) + 
    facet_grid(Be+Bu~Ne+Nu, labeller = label_parsed) + 
    theme_bw() + 
    xlab("Power") + 
    ylab("Empirical FDR") +
    theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        # axis.text.x = element_blank(),
        #  axis.ticks = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        title = element_text(size = 14),
        legend.position = "top") +
    scale_color_manual(name = "", labels=c("GEM", "Clump", "COJO", "SharePro(variant)", "SharePro(effect)"),
                                values=c("gold", "chartreuse2", "darkgreen", "purple", "darkblue")) +
    scale_shape_manual(name = "", labels=c("GEM", "Clump", "COJO", "SharePro(variant)", "SharePro(effect)"),
                                values=c(15,16,17,9,10)) +
    scale_y_continuous(breaks = c(0,0.5,1)) + 
    scale_x_continuous(breaks = c(0,0.5,1)) -> plt
ggsave(paste0(fdir,"sim_gxe_3.pdf"),plt,height = 10, width = 10)

ggplot(df_pwr_fdr[df_pwr_fdr$K==2,], aes(x=Power, y=FDR, color=Method, shape=Method)) + 
    geom_point(size=3) + 
    facet_grid(Be+Bu~Ne+Nu, labeller = label_parsed) + 
    theme_bw() + 
    xlab("Power") + 
    ylab("Empirical FDR") +
    theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        # axis.text.x = element_blank(),
        #  axis.ticks = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        title = element_text(size = 14),
        legend.position = "top") +
    scale_color_manual(name = "", labels=c("GEM", "Clump", "COJO", "SharePro(variant)", "SharePro(effect)"),
                                values=c("gold", "chartreuse2", "darkgreen", "purple", "darkblue")) +
    scale_shape_manual(name = "", labels=c("GEM", "Clump", "COJO", "SharePro(variant)", "SharePro(effect)"),
                                values=c(15,16,17,9,10)) +
    scale_y_continuous(breaks = c(0,0.5,1)) + 
    scale_x_continuous(breaks = c(0,0.5,1)) -> plt
ggsave(paste0(fdir,"sim_gxe_2.pdf"),plt,height = 10, width = 10)

ggplot(df_pwr_fdr[df_pwr_fdr$K==1,], aes(x=Power, y=FDR, color=Method, shape=Method)) + 
    geom_point(size=3) + 
    facet_grid(Be+Bu~Ne+Nu, labeller = label_parsed) + 
    theme_bw() + 
    xlab("Power") + 
    ylab("Empirical FDR") +
    theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        # axis.text.x = element_blank(),
        #  axis.ticks = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        title = element_text(size = 14),
        legend.position = "top") +
    scale_color_manual(name = "", labels=c("GEM", "Clump", "COJO", "SharePro(variant)", "SharePro(effect)"),
                                values=c("gold", "chartreuse2", "darkgreen", "purple", "darkblue")) +
    scale_shape_manual(name = "", labels=c("GEM", "Clump", "COJO", "SharePro(variant)", "SharePro(effect)"),
                                values=c(15,16,17,9,10)) +
    scale_y_continuous(breaks = c(0,0.5,1)) + 
    scale_x_continuous(breaks = c(0,0.5,1)) -> plt
ggsave(paste0(fdir,"sim_gxe_1.pdf"),plt,height = 10, width = 10)

load("../doc/calibration.Rdata")
stats = data.frame(setting = rep(names(nsharepro_pval_lst),each = 3),significance = c(0.1,0.05,0.01),type_1_error_rate_GEM = NA,type_1_error_rate_SharePro = NA,type_1_error_rate_COJO = NA,type_1_error_rate_Clump = NA)
for (j in 1:nrow(stats)) {
    name = stats$setting[j]
    thres = stats$significance[j]
    ps = ngem_pval_lst[[name]]
    stats$type_1_error_rate_GEM[j] = sum(ps < thres) / length(ps)
    ps = ngem_pval_lst[[name]][nsharepro_pval_lst[[name]]==1]
    stats$type_1_error_rate_SharePro[j] = sum(ps < thres) / length(ps)
    ps = ngem_pval_lst[[name]][ncojo_pval_lst[[name]]==1]
    stats$type_1_error_rate_COJO[j] = sum(ps < thres) / length(ps)
    ps = ngem_pval_lst[[name]][nclump_pval_lst[[name]]==1]
    stats$type_1_error_rate_Clump[j] = sum(ps < thres) / length(ps)
}

stats_formatted = data.frame(Number_of_causal_variants=unlist(lapply(strsplit(stats$setting, '_'), function(x){x[1]})),
                            Sample_size_in_exposed=25000,
                            Sample_size_in_unexposed=rep(c(25000,10000,5000), each=3),
                            Effect_size_in_exposed=0.05,
                            Effect_size_in_unexposed=0.05,
                            Significance_threshold=stats$significance,
                            Type_1_error_rate_GEM=stats$type_1_error_rate_GEM,
                            Type_1_error_rate_Clump=stats$type_1_error_rate_Clump,
                            Type_1_error_rate_COJO=stats$type_1_error_rate_COJO,
                            Type_1_error_rate_SharePro=stats$type_1_error_rate_SharePro)

write.table(stats_formatted, file = '../doc/calibration_formatted.txt', quote=F, sep='\t', row.names=F, col.names=T)