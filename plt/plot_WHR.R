library(ggplot2)
library(ggrepel)
set.seed(2023)

cs <- read.table("../doc/WHR_sex_hormone_LDL_TG_HDL.txt",sep = "\t",header = T)

cs$significance <- ifelse(cs$p_diff < 5e-8, "Genome-wide GxSex in WHRadjBMI",
                          ifelse(cs$p_diff < 0.05 / nrow(cs), "Effect group-level GxSex in WHRadjBMI", "Non-significant"))
cs$significance <- factor(cs$significance, levels = c("Genome-wide GxSex in WHRadjBMI",
                                                      "Effect group-level GxSex in WHRadjBMI",
                                                      "Non-significant"))
cs$p_diff_Testo <- -log10(pchisq(((cs$bF_Testo - cs$bM_Testo) / sqrt(cs$seF_Testo^2 + cs$seM_Testo^2))^2, df = 1, lower.tail = F))
cs$p_diff_SHBG <- -log10(pchisq(((cs$bF_SHBG - cs$bM_SHBG) / sqrt(cs$seF_SHBG^2 + cs$seM_SHBG^2))^2, df = 1, lower.tail = F))
cs$p_diff_LDL <- -log10(pchisq(((cs$bF_LDL - cs$bM_LDL) / sqrt(cs$seF_LDL^2 + cs$seM_LDL^2))^2, df = 1, lower.tail = F))
cs$p_diff_TG <- -log10(pchisq(((cs$bF_TG - cs$bM_TG) / sqrt(cs$seF_TG^2 + cs$seM_TG^2))^2, df = 1, lower.tail = F))
cs$p_diff_HDL <- -log10(pchisq(((cs$bF_HDL - cs$bM_HDL) / sqrt(cs$seF_HDL^2 + cs$seM_HDL^2))^2, df = 1, lower.tail = F))

cs$p_diff_Testo[cs$p_diff_Testo > 10] <- 10
cs$p_diff_SHBG[cs$p_diff_SHBG > 10] <- 10
cs$p_diff_LDL[cs$p_diff_LDL > 10] <- 10
cs$p_diff_TG[cs$p_diff_TG > 10] <- 10
cs$p_diff_HDL[cs$p_diff_HDL > 10] <- 10


cs$diff_beta_Testo <- abs(cs$bF_Testo * sqrt(2 * cs$mfF_Testo * (1 - cs$mfF_Testo)) - cs$bM_Testo * sqrt(2 * cs$mfM_Testo * (1 - cs$mfM_Testo)))
cs$diff_beta_SHBG <- abs(cs$bF_SHBG * sqrt(2 * cs$mfF_SHBG * (1 - cs$mfF_SHBG)) - cs$bM_SHBG * sqrt(2 * cs$mfM_SHBG * (1 - cs$mfM_SHBG)))
cs$diff_beta_LDL <-abs(cs$bF_LDL * sqrt(2 * cs$mfF_LDL * (1 - cs$mfF_LDL)) - cs$bM_LDL * sqrt(2 * cs$mfM_LDL * (1 - cs$mfM_LDL)))
cs$diff_beta_HDL <-abs(cs$bF_HDL * sqrt(2 * cs$mfF_HDL * (1 - cs$mfF_LDL)) - cs$bM_HDL * sqrt(2 * cs$mfM_HDL * (1 - cs$mfM_HDL)))
cs$diff_beta_TG <- abs(cs$bF_TG * sqrt(2 * cs$mfF_TG * (1 - cs$mfF_TG)) - cs$bM_TG * sqrt(2 * cs$mfM_TG * (1 - cs$mfM_TG)))
cs$diff_beta_WHRadjBMI <- abs(cs$betaF-cs$betaM)

cs$label <- ifelse(cs$p_diff_Testo > -log10(0.05/nrow(cs)),cs$Nearest_gene,"")
ggplot(cs, aes(x = diff_beta_Testo, 
               y = diff_beta_WHRadjBMI, 
               color = significance,
               size = p_diff_Testo,
               label = label)) +
  geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  xlab(expression(Testosterone:"|"~beta[F]-beta[M]~"|")) +
  ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
  labs(color = "",
       size = expression(-log[10](p[GxSex~Testosterone]))) +
  scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
  scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
  geom_text_repel(show.legend = F) +
  guides(size = "legend",color = "none")-> plt
ggsave("../Fig/WHR_Testo.pdf",plt,height = 6,width = 6.5)

cs$label <- ifelse(cs$p_diff_SHBG > -log10(0.05/nrow(cs)),cs$Nearest_gene,"")
ggplot(cs, aes(x = diff_beta_SHBG, 
               y = diff_beta_WHRadjBMI, 
               color = significance,
               size = p_diff_SHBG,
               label = label)) +
  geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  xlab(expression(SHBG:"|"~beta[F]-beta[M]~"|")) +
  ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
  labs(color = "",
       size = expression(-log[10](p[GxSex~SHBG]))) +
  scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
  scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
  geom_text_repel(show.legend = F) +
  guides(size = "legend",color = "none")-> plt
ggsave("../Fig/WHR_SHBG.pdf",plt,height = 6,width = 6.5)

cs$label <- ifelse(cs$p_diff_LDL > -log10(0.05/nrow(cs)),cs$Nearest_gene,"")
ggplot(cs, aes(x = diff_beta_LDL, 
               y = diff_beta_WHRadjBMI, 
               color = significance,
               size = p_diff_LDL,
               label = label)) +
  geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  xlab(expression(LDL:"|"~beta[F]-beta[M]~"|")) +
  ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
  labs(color = "",
       size = expression(-log[10](p[GxSex~LDL]))) +
  scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
  scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
  geom_text_repel(show.legend = F) +
  guides(size = "legend",color = "none")-> plt
ggsave("../Fig/WHR_LDL.pdf",plt,height = 6,width = 6.5)


cs$label <- ifelse(cs$p_diff_HDL > -log10(0.05/nrow(cs)),cs$Nearest_gene,"")
ggplot(cs, aes(x = diff_beta_HDL, 
               y = diff_beta_WHRadjBMI, 
               color = significance,
               size = p_diff_HDL,
               label = label)) +
  geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  xlab(expression(HDL:"|"~beta[F]-beta[M]~"|")) +
  ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
  labs(color = "",
       size = expression(-log[10](p[GxSex~HDL]))) +
  scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
  scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
  geom_text_repel(show.legend = F) +
  guides(size = "legend",color = "none")-> plt
ggsave("../Fig/WHR_HDL.pdf",plt,height = 6,width = 6.5)

cs$label <- ifelse(cs$p_diff_TG > -log10(0.05/nrow(cs)),cs$Nearest_gene,"")
ggplot(cs, aes(x = diff_beta_TG, 
               y = diff_beta_WHRadjBMI, 
               color = significance,
               size = p_diff_TG,
               label = label)) +
  geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  xlab(expression(TG:"|"~beta[F]-beta[M]~"|")) +
  ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
  labs(color = "",
       size = expression(-log[10](p[GxSex~TG]))) +
  scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
  scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
  geom_text_repel(show.legend = F) +
  guides(size = "legend",color = "none")-> plt
ggsave("../Fig/WHR_TG.pdf",plt,height = 6,width = 6.5)

cs$label <- ifelse(cs$pAFdiff > -log10(0.05/nrow(cs)),cs$Nearest_gene,"")
ggplot(cs, aes(x = AFdiff, 
               y = diff_beta_WHRadjBMI, 
               color = significance,
               size = pAFdiff,
               label = label)) +
  geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  xlab(expression("|"~MAF[F]-MAF[M]~"|")) +
  ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
  labs(color = "",
       size = expression(-log[10](p[sex~diff.~MAF]))) +
  scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
  scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
  xlim(0,0.005) +
  geom_text_repel(show.legend = F) +
  guides(size = "legend",color = "none")-> plt
ggsave("../Fig/WHR_MAF.pdf",plt,height = 6,width = 6.5)

ggplot(cs, aes(x = AFdiff, 
               y = diff_beta_WHRadjBMI, 
               color = significance,
               size = pAFdiff,
               label = label)) +
  geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  xlab(expression("|"~MAF[F]-MAF[M]~"|")) +
  ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
  labs(color = "",
       size = expression(-log[10](p[sex~diff.~MAF]))) +
  scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
  scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
  xlim(0,0.005) +
  geom_text_repel(show.legend = F) +
  guides(color = "legend", size="none")-> plt
ggsave("../Fig/WHR_legend.pdf",plt,height = 7,width = 16)

## Correlation test

cor.test(cs$diff_beta_WHRadjBMI, cs$diff_beta_SHBG, method = "spearman")
cor.test(cs$diff_beta_WHRadjBMI, cs$diff_beta_Testo, method = "spearman")
cor.test(cs$diff_beta_WHRadjBMI, cs$diff_beta_TG, method = "spearman")
cor.test(cs$diff_beta_WHRadjBMI, cs$diff_beta_LDL, method = "spearman")
cor.test(cs$diff_beta_WHRadjBMI, cs$diff_beta_HDL, method = "spearman")
cor.test(cs$diff_beta_WHRadjBMI, cs$AFdiff, method = "spearman")

