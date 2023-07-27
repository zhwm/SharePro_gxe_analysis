library(PRROC)
library(ggplot2)
library(data.table)

wdir = "../sim/"
odir = "../doc/"
fdir = "../Fig/"
setwd(wdir)

### An example ###############
# Locus1/21
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
SH = fread(paste0(wdir,lo,folder,'SH/',substr(rt,1,1),x,'.z_',substr(rt,2,2),x,'.z.cs'))
SHp = fread(paste0(wdir,lo,folder,'SH/',substr(rt,1,1),x,'.z_',substr(rt,2,2),x,'.z.snp'))
SPA = fread(paste0(wdir,lo,folder,'SA/A',rt,x,'.z.cs'))
SPAp = fread(paste0(wdir,lo,folder,'SA/A',rt,x,'.z.snp'))
gem$shp <- NA
cslist = strsplit(SH$cs,'/')
for (i in 1:nrow(SH)) {
  e = cslist[i]
  gem$shp[gem$SNPID %in% unlist(e)]=SH$p_diff[i]
}

pvalplt = data.frame(chr = rep(gem$CHR,9),
                     pos = rep(gem$POS,9),
                     P = c(plink$pval1,
                           plink$pval2,
                           gem$P_Value_Interaction,
                           p.adjust(gem$P_Value_Interaction,'bonferroni'),
                           p.adjust(gem$P_Value_Interaction,'BH'),
                           plink$P_GXE,
                           gem$P_Value_Marginal,
                           gem$P_Value_Joint,
                           gem$shp),
                     LD = rep(ld[,which(gem$SNPID==csnp[x,])]^2,9),
                     Adjustment = rep(c("Exposed association test",
                                        "Unexposed association test",
                                        "GEM interaction test",
                                        "GEM interaction test (Bonferroni-adjusted)",
                                        "GEM interaction test (BH-adjusted)",
                                        "PLINK interaction test",
                                        "Combined 1-df association test",
                                        "Combined 2-df association test",
                                        "SharePro interaction test"),each=nrow(gem)))
## association test
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
ggsave(paste0(fdir, "E_UE_association_example.pdf"),plt, height = 5, width = 4)

#GEM GxE test
ggplot(pvalplt[pvalplt$Adjustment=="GEM interaction test" | pvalplt$Adjustment=="GEM interaction test (BH-adjusted)",], aes(x = pos / 1000000, y = -log10(P), color=LD)) +
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
ggsave(paste0(fdir, "GEM_GxE_example.pdf"),plt,height = 5, width = 4)

#SharePro GxE test
ggplot(pvalplt[pvalplt$Adjustment=="SharePro interaction test",], aes(x = pos / 1000000, y = -log10(P), color=LD)) +
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
ggsave(paste0(fdir, "SharePro_GxE_example.pdf"),plt,height = 3, width = 4)

pipplt = data.frame(chr = rep(gem$CHR,4),
                    pos = rep(gem$POS,4),
                    P = c(SPAp$vProb,
                          SHp$vProb),
                    LD = rep(ld[,which(gem$SNPID==csnp[x,])]^2,4),
                    Adjustment = rep(c("Combined fine-mapping","SharePro fine-mapping"),each=nrow(gem)))

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
ggsave(paste0(fdir, "SharePro_finemap_example.pdf"),plt,height = 3, width = 4)

#################################################

## Summary 
get_cs_power_fdr <- function(cs,tsnp){
  if (length(cs)==0) {
    power = 0
    fdr = 0
  } else {
    if (length(tsnp)==0) {
      power = 0
      fdr = 1
    } else {
      tp = length(intersect(cs,tsnp))
      power = tp / length(tsnp)
      fdr = 1 - tp / length(cs)
    }
  }
  return(c(power,fdr))
}

get_auprc_vec <- function(pip,gt){
  cv = pr.curve(pip[gt],pip[!gt])
  return(cv$auc.integral)
}

ctf <- 0.01

auprc_vec <- c() # AUPRC in identifying causal variants
pval_pr_vec <- c() # Precision recall on GxE

folderlist <- c('1_0.05_-0.05/','1_0.05_-0.02/','1_0.05_0/','1_0.05_0.02/','1_0.05_0.05/',
                '2_0.05_-0.05/','2_0.05_-0.02/','2_0.05_0/','2_0.05_0.02/','2_0.05_0.05/',
                '3_0.05_-0.05/','3_0.05_-0.02/','3_0.05_0/','3_0.05_0.02/','3_0.05_0.05/')
Locilist <- c("Locus1/","Locus2/","Locus3/")

for (folder in folderlist) {
  params = unlist(strsplit(gsub('/','',folder),'_'))
  for (rt in c('CL','CM','CS')) {
    for (lo in Locilist) {
      gem_pmv <- c()
      gem_pjv <- c()
      spav <- c()
      shv <- c()
      GTv <- c()
      csnp = read.table(paste0(wdir,lo,folder,'csnp.txt'),header=F)
      print(paste0(lo,folder,rt))
      for (x in 1:50) {
        gem = fread(paste0(wdir,lo,folder,rt,x))
        plink = fread(paste0(wdir,lo,folder,rt,x,'.qassoc.gxe'))
        SH = fread(paste0(wdir,lo,folder,'SH/',substr(rt,1,1),x,'.z_',substr(rt,2,2),x,'.z.cs'))
        SPA = fread(paste0(wdir,lo,folder,'SA/A',rt,x,'.z.cs'))
        SHp = fread(paste0(wdir,lo,folder,'SH/',substr(rt,1,1),x,'.z_',substr(rt,2,2),x,'.z.snp'))
        SPAp = fread(paste0(wdir,lo,folder,'SA/A',rt,x,'.z.snp'))
        spav <- c(spav,SPAp$vProb)
        shv <- c(shv,SHp$vProb)
        GTv <- c(GTv,gem$SNPID %in% csnp[x,])
        gem_pmv <- c(gem_pmv,-gem$P_Value_Marginal)
        gem_pjv <- c(gem_pjv,-gem$P_Value_Joint)
        
        if (params[2]==params[3]) {
          print("No gxe variants")
          gxesnp = c()
        } else {
          gxesnp = csnp[x,]
        }
        
        gem$shp <- NA
        totalp = 0
        SHdiff = SH[SH$p_diff<ctf,]
        if (nrow(SH)>0) {
          cslist = strsplit(SH$cs,'/')
          for (i in 1:nrow(SH)) {
            e = cslist[i]
            gem$shp[gem$SNPID %in% unlist(e)]=SH$p_diff[i]
            if (length(intersect(unlist(e),gxesnp))==1 & SH$p_diff[i]<ctf) {
              totalp = totalp + 1
            }
          }
        }
        
        if (nrow(SHdiff)==0) {
          powercs = 0
          fdrcs = 0
        } else {
          if (length(gxesnp)==0) {
            powercs = 0
            fdrcs = 1
          } else {
            powercs = totalp / length(gxesnp)
            fdrcs = 1 - totalp / nrow(SHdiff)
          }
        }
        
        pval_pr_vec <- c(pval_pr_vec,c(powercs,fdrcs))
        pvalues <- list(gem$P_Value_Interaction,
                        p.adjust(gem$P_Value_Interaction,'bonferroni'),
                        p.adjust(gem$P_Value_Interaction,'BH'),
                        plink$P_GXE,
                        gem$shp)
        pval_pr_vec <- c(pval_pr_vec, unlist(lapply(pvalues, function(y)
        {get_cs_power_fdr(gem$SNPID[y<ctf & !is.na(y)],gxesnp)})))
      }
      pips <- list(spav,
                   shv,
                   gem_pmv,
                   gem_pjv)
      auprc_vec <- c(auprc_vec,unlist(lapply(pips, function(y){get_auprc_vec(y,GTv)})))
    }
  }
}



######## Record data
df_auprc <- data.frame(auprc = auprc_vec,
                       method = c("Combined fine-mapping",
                                  "SharePro fine-mapping",
                                  "1-df association",
                                  "2-df association"),
                       rt = rep(rep(c("CL",
                                      "CM",
                                      "CS"),each=4*length(Locilist)),length(folderlist)),
                       folder = rep(folderlist,each=4*3*length(Locilist)))
df_auprc$K  = unlist(lapply(strsplit(df_auprc$folder,'_'),function(x){x[1]}))
df_auprc$Ne = 25000
df_auprc$Nu = factor(unlist(lapply(df_auprc$rt, function(x){strsplit(x,2,2)})),
                     labels = c('25000','10000','5000'))
df_auprc$Be = unlist(lapply(strsplit(df_auprc$folder,'_'),function(x){x[2]}))
df_auprc$Bu = unlist(lapply(strsplit(df_auprc$folder,'_'),function(x){gsub('/','',x[3])}))

write.table(df_auprc[, c('K', 'Ne', 'Nu', 'Be', 'Bu', 'method', 'auprc')], paste0(odir, "SharePro_gxe_sim_auprc.csv"),sep=',',row.names=F,col.names=T,quote=F)

df_auprc_median = aggregate(auprc~K+Ne+Nu+Be+Bu+method,data = df_auprc, median)
df_auprc_median = reshape(df_auprc_median, direction = 'wide', idvar=c('K', 'Ne', 'Nu', 'Be', 'Bu'), timevar='method')
colnames(df_auprc_median) = c('K', 'Ne', 'Nu', 'Be', 'Bu', '1-df association', '2-df association', 'Combined fine-mapping', 'SharePro fine-mapping')
write.table(df_auprc_median, paste0(odir, "SharePro_gxe_sim_auprc_median.csv"), sep=',',row.names=F,col.names=T,quote=F)

df_gxep_pr <- data.frame(power = pval_pr_vec[seq(1, length(pval_pr_vec), 2)],
                         fdr = pval_pr_vec[seq(2, length(pval_pr_vec), 2)],
                         method = c("SharePro (effect)",
                                    "GEM",
                                    "GEM-Bonferroni",
                                    "GEM-BH",
                                    "PLINK",
                                    "SharePro (variant)"),
                         rt = rep(rep(c("CL",
                                        "CM",
                                        "CS"),each=50*length(Locilist)*6),length(folderlist)),
                         folder = rep(folderlist,each=3*6*50*length(Locilist)))

df_gxep_pr$K  = unlist(lapply(strsplit(df_gxep_pr$folder,'_'),function(x){x[1]}))
df_gxep_pr$Ne = 25000
df_gxep_pr$Nu = factor(unlist(lapply(df_gxep_pr$rt, function(x){strsplit(x,2,2)})),
                       labels = c('25000','10000','5000'))
df_gxep_pr$Be = unlist(lapply(strsplit(df_gxep_pr$folder,'_'),function(x){x[2]}))
df_gxep_pr$Bu = unlist(lapply(strsplit(df_gxep_pr$folder,'_'),function(x){gsub('/','',x[3])}))

write.table(df_gxep_pr[, c('K', 'Ne', 'Nu', 'Be', 'Bu', 'method', 'power', 'fdr')], paste0(odir, "SharePro_gxe_sim_gxep_pr.csv"),sep=',',row.names=F,col.names=T,quote=F)



######## AUPRC plot
df_auprc = read.csv(paste0(odir, "SharePro_gxe_sim_auprc.csv"), header=T)

df_auprc$K = paste0("K:", df_auprc$K)
df_auprc$Ne = paste0("N[e]:", df_auprc$Ne)
df_auprc$Nu = paste0("N[u]:", df_auprc$Nu)
df_auprc$Nu = factor(df_auprc$Nu, levels = c("N[u]:25000", "N[u]:10000", "N[u]:5000"))
df_auprc$Be = paste0("beta[e]:", df_auprc$Be)
df_auprc$Bu = paste0("beta[u]:", df_auprc$Bu)
df_auprc$Bu = factor(df_auprc$Bu, levels = c("beta[u]:-0.05", "beta[u]:-0.02", "beta[u]:0", "beta[u]:0.02", "beta[u]:0.05"))


ggplot(df_auprc[df_auprc$K=="K:3",], aes(x=method, y=auprc, color=method)) + 
  facet_grid(K+Be+Bu~Nu+Ne,labeller = label_parsed) +
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
  scale_color_manual(values = c("orange","brown","royalblue","darkblue")) + 
  labs(color = "") -> plt
ggsave(paste0(fdir, "SharePro_gxe_sim_auprc.pdf"),plt,height = 15, width = 10)


####### GxE plot
df_gxep_pr = read.csv(paste0(odir, "SharePro_gxe_sim_gxep_pr.csv"), header=T)
df_gxe_summary = aggregate(cbind(power, fdr) ~ K+Ne+Nu+Be+Bu+method, data = df_gxep_pr, mean)
write.table(df_gxe_summary,paste0(odir, "SharePro_gxe_sim_gxe_summary.csv"),sep=',',row.names=F,col.names=T,quote=F)
df_gxe_summary$K = paste0("K:", df_gxe_summary$K)
df_gxe_summary$Ne = paste0("N[e]:", df_gxe_summary$Ne)
df_gxe_summary$Nu = paste0("N[u]:", df_gxe_summary$Nu)
df_gxe_summary$Nu = factor(df_gxe_summary$Nu, levels = c("N[u]:25000", "N[u]:10000", "N[u]:5000"))
df_gxe_summary$Be = paste0("beta[e]:", df_gxe_summary$Be)
df_gxe_summary$Bu = paste0("beta[u]:", df_gxe_summary$Bu)
df_gxe_summary$Bu = factor(df_gxe_summary$Bu, levels = c("beta[u]:-0.05", "beta[u]:-0.02", "beta[u]:0", "beta[u]:0.02", "beta[u]:0.05"))
df_gxe_summary$method = factor(df_gxe_summary$method,
                               c("GEM",
                                 "PLINK",
                                 "GEM-Bonferroni",
                                 "GEM-BH",
                                 "SharePro (variant)",
                                 "SharePro (effect)"))


ggplot(df_gxe_summary[df_gxe_summary$K==paste0('K:',3),],
       aes(x=power,y=1-fdr,color=method,shape=method)) + 
  geom_point(size=3) +
  facet_grid(K+Be+Bu~Nu+Ne,labeller = label_parsed) +
  theme_bw() +
  xlab("Power") + 
  ylab("1-FDR") +
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
  scale_color_manual(name = "", labels=c("GEM",
                                         "PLINK",
                                         "GEM-Bonferroni",
                                         "GEM-BH",
                                         "SharePro (variant)",
                                         "SharePro (effect)"),values = c("gold","darkgoldenrod1","chartreuse2","darkgreen","purple","darkblue")) + 
  scale_shape_manual(name = "",labels=c("GEM",
                                        "PLINK",
                                        "GEM-Bonferroni",
                                        "GEM-BH",
                                        "SharePro (variant)",
                                        "SharePro (effect)"),values = c(15:18,9,10)) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  scale_x_continuous(breaks = c(0,0.5,1)) -> plt
ggsave(paste0(fdir,"SharePro_gxe_sim_gxe_summary.pdf"),plt,height = 15, width = 10)

ggplot(df_gxe_summary[df_gxe_summary$K==paste0('K:',3),],
       aes(x=method,y=fdr,color=method, shape=method)) + 
  geom_point(size = 3) +
  facet_grid(K+Be+Bu~Nu+Ne,labeller = label_parsed) +
  theme_bw() + 
  xlab("Method") + 
  ylab("FDR") +
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
  scale_color_manual(name = "", labels=c("GEM",
                                         "PLINK",
                                         "GEM-Bonferroni",
                                         "GEM-BH",
                                         "SharePro (variant)",
                                         "SharePro (effect)"),values = c("gold","darkgoldenrod1","chartreuse2","darkgreen","purple","darkblue")) + 
  scale_shape_manual(name = "",labels=c("GEM",
                                        "PLINK",
                                        "GEM-Bonferroni",
                                        "GEM-BH",
                                        "SharePro (variant)",
                                        "SharePro (effect)"),values = c(15:18,9,10)) +
  scale_y_continuous(breaks = c(0,0.5,1)) -> plt
ggsave(paste0(fdir,"SharePro_gxe_sim_gxe_summary_fdr.pdf"),plt,height = 15, width = 10)

#####################################################
## time 
timedt = read.csv(paste0(odir,'time.summary'),header=F,sep=' ')
timedt$method = unlist(lapply(strsplit(timedt$V1,'/'), function(x){x[5]}))
timedf = timedt[timedt$method=="gem.time" | timedt$method=="plink.time" | timedt$method=="CL.time",]
timedf$method[timedf$method=="CL.time"]="SharePro"
timedf$method[timedf$method=="plink.time"]="PLINK"
timedf$method[timedf$method=="gem.time"]="GEM"
timedf$second <- timedf$V2/50

dftime = aggregate(second~method,data=timedf,mean)
dftime$sd = aggregate(second~method,data=timedf,sd)$second
colnames(dftime) <- c("Method","Run time average(s)", "Run time standard deviation(s)")
write.table(dftime,paste0(odir, 'time.csv'),sep=',',col.names = T,row.names = F)

