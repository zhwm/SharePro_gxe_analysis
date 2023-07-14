library(ggplot2)
library(ramwas)
library(data.table)

gxesnp <- fread('../dat/GxE/FFR_smoke_gxesnp.txt')

pdf(file = "../Fig/FFR_smoking_CN_gxe.pdf", height = 4, width = 6)
cnplt <- manPlotPrepare(gxesnp$pCN,
                        gxesnp$CHR,
                        gxesnp$POS)
manPlotFast(cnplt,ylim = c(0,18),colorSet = c("plum","darkblue"))
abline(h = -log10(5e-8),lty = 2,col = "red")
dev.off()

pdf(file = "../Fig/FFR_smoking_CE_gxe.pdf", height = 4, width = 6)
ceplt <- manPlotPrepare(gxesnp$pCE,
                        gxesnp$CHR,
                        gxesnp$POS)
manPlotFast(ceplt,ylim = c(0,18),colorSet = c("plum","darkblue"))
abline(h = -log10(5e-8),lty = 2,col = "red")
dev.off()

pdf(file = "../Fig/FFR_smoking_EN_gxe.pdf", height = 4, width = 6)
enplt <- manPlotPrepare(gxesnp$pEN,
                        gxesnp$CHR,
                        gxesnp$POS)
manPlotFast(enplt,ylim = c(0,18),colorSet = c("plum","darkblue"))
abline(h = -log10(5e-8),lty = 2,col = "red")
dev.off()

cs <- read.table("../doc/FFR_smoke_CEN_cgene.txt",sep = "\t",header = T)
cen_n <- c(27409,
           100285,
           155951)
cs$betaC <- as.numeric(unlist(lapply(strsplit(cs$beta,','), function(x){x[1]})))
cs$betaE <- as.numeric(unlist(lapply(strsplit(cs$beta,','), function(x){x[2]})))
cs$betaN <- as.numeric(unlist(lapply(strsplit(cs$beta,','), function(x){x[3]})))

cs$pCE <- pchisq((cs$betaC - cs$betaE)**2 / (1/cen_n[1] + 1/cen_n[2]), 1, lower.tail = F)
cs$pEN <- pchisq((cs$betaE - cs$betaN)**2 / (1/cen_n[2] + 1/cen_n[3]), 1, lower.tail = F)
cs$pCN <- cs$p_diff

write.table(cs[,c('Top_variant', 'betaC', 'betaE', 'betaN', 'pCE', 'pCN', 'pEN', 'Nearest_gene', 
                'Distance_to_nearest_gene','eQTL', 'sQTL', 'cs', 'rsid')], 
          file='../doc/FFR_smoking_gxe_summary.txt', sep='\t', row.names = F)

pctf = 0.05/nrow(cs)/3

gxe = merge(y=gxesnp[,c('SNP','CHR','POS')], x=cs[,c('Top_variant','pCE','pCN','pEN')], by.y="SNP", by.x="Top_variant")

pdf(file = "../Fig/FFR_smoking_CN_SharePro.pdf", height = 4, width = 6)
cn <- manPlotPrepare(gxe$pCN,
                     as.factor(gxe$CHR),
                     gxe$POS)
manPlotFast(cn,ylim = c(0,18),colorSet = c("plum","darkblue"))
abline(h = -log10(pctf),lty = 2,col = "red")
dev.off()

pdf(file = "../Fig/FFR_smoking_CE_SharePro.pdf", height = 4, width = 6)
ce <- manPlotPrepare(gxe$pCE,
                     as.factor(gxe$CHR),
                     gxe$POS)
manPlotFast(ce,ylim = c(0,18),colorSet = c("plum","darkblue"))
abline(h = -log10(pctf),lty = 2,col = "red")
dev.off()

pdf(file = "../Fig/FFR_smoking_EN_SharePro.pdf", height = 4, width = 6)
en <- manPlotPrepare(gxe$pEN,
                     as.factor(gxe$CHR),
                     gxe$POS)
manPlotFast(en,ylim = c(0,18),colorSet = c("plum","darkblue"))
abline(h = -log10(pctf),lty = 2,col = "red")
dev.off()

sigcs = cs[cs$pCE<pctf | cs$pCN<pctf | cs$pEN<pctf,]
pltdat <- data.frame(Gene = sigcs$Nearest_gene,
                     b = c(sigcs$betaC,sigcs$betaE, sigcs$betaN),
                     SE = rep(sqrt(1/cen_n), each=nrow(sigcs)),
                     Population = rep(c("Current smokers", "Past smokers", "Never smokers"),each = nrow(sigcs)))
pltdat$Population <- factor(pltdat$Population, levels = c("Current smokers", "Past smokers", "Never smokers"))
pltdat$Gene <- factor(pltdat$Gene,levels = rev(c('CHRNA3', 'ADAM19', 'UBR1')))
ggplot(pltdat, aes(x = b, y = Gene, color = Population)) +
  geom_vline(xintercept = 0, lty = 2, color = "darkgrey") +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_errorbarh(aes(xmin = b - 1.96 * SE,
                     xmax = b + 1.96 * SE),height = 0,position = position_dodge(width = 0.3)) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.position = "right") +
  scale_color_manual(values = c("darkred", "orange", "darkgreen")) +
  labs(color = "") +
  xlab(expression(beta[FFR])) +
  ylab("Nearest gene") -> plt
ggsave("../Fig/FFR_smoking_gene.pdf", plt, height = 6, width = 6)
