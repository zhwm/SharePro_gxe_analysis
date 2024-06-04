library(data.table)
library(ggplot2)
library(ggrepel)

WHRadjBMI_cs = fread('../doc/WHR_summary.txt')
WHRadjBMI_cs$significance = factor(WHRadjBMI_cs$significance, levels = c("Genome-wide GxSex in WHRadjBMI","Effect group-level GxSex in WHRadjBMI","Non-significant"))

cor.test(WHRadjBMI_cs$diff_beta_WHRadjBMI, WHRadjBMI_cs$diff_beta_SHBG, method = "spearman")
cor.test(WHRadjBMI_cs$diff_beta_WHRadjBMI, WHRadjBMI_cs$diff_beta_T, method = "spearman")
cor.test(WHRadjBMI_cs$diff_beta_WHRadjBMI, WHRadjBMI_cs$diff_beta_LDL, method = "spearman")
cor.test(WHRadjBMI_cs$diff_beta_WHRadjBMI, WHRadjBMI_cs$diff_beta_TG, method = "spearman")
cor.test(WHRadjBMI_cs$diff_beta_WHRadjBMI, WHRadjBMI_cs$diff_beta_HDL, method = "spearman")
cor.test(WHRadjBMI_cs$diff_beta_WHRadjBMI, WHRadjBMI_cs$diff_maf, method = "spearman")

WHRadjBMI_cs$label = ifelse(WHRadjBMI_cs$p_diff_SHBG > -log10(0.05/nrow(WHRadjBMI_cs)), WHRadjBMI_cs$NearestGene, "")
ggplot(WHRadjBMI_cs, aes(x = diff_beta_SHBG, y = diff_beta_WHRadjBMI, color = significance, size = p_diff_SHBG, label = label)) +
    geom_point() +
    theme_bw() +
    theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
    xlab(expression(SHBG:"|"~beta[F]-beta[M]~"|")) +
    ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
    labs(color = "", size = expression(-log[10](p[GxSex~SHBG]))) +
    scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
    scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
    geom_text_repel(show.legend = F, max.overlaps = 50) +
    guides(size = "legend",color = "none") +
    ylim(0,0.004) -> plt
ggsave("../fig/WHR_SHBG.pdf",plt,height = 6,width = 6.5)

WHRadjBMI_cs$label = ifelse(WHRadjBMI_cs$p_diff_T > -log10(0.05/nrow(WHRadjBMI_cs)), WHRadjBMI_cs$NearestGene, "")
ggplot(WHRadjBMI_cs, aes(x = diff_beta_T, y = diff_beta_WHRadjBMI, color = significance, size = p_diff_T, label = label)) +
    geom_point() +
    theme_bw() +
    theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
    xlab(expression(Testosterone:"|"~beta[F]-beta[M]~"|")) +
    ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
    labs(color = "", size = expression(-log[10](p[GxSex~Testosterone]))) +
    scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
    scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
    geom_text_repel(show.legend = F, max.overlaps = 50) +
    guides(size = "legend",color = "none") +
    ylim(0,0.004) -> plt
ggsave("../fig/WHR_Testo.pdf",plt,height = 6,width = 6.5)

WHRadjBMI_cs$label = ifelse(WHRadjBMI_cs$p_diff_LDL > -log10(0.05/nrow(WHRadjBMI_cs)), WHRadjBMI_cs$NearestGene, "")
ggplot(WHRadjBMI_cs, aes(x = diff_beta_LDL, y = diff_beta_WHRadjBMI, color = significance, size = p_diff_LDL, label = label)) +
    geom_point() +
    theme_bw() +
    theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
    xlab(expression(LDL:"|"~beta[F]-beta[M]~"|")) +
    ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
    labs(color = "", size = expression(-log[10](p[GxSex~LDL]))) +
    scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
    scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
    geom_text_repel(show.legend = F, max.overlaps = 50) +
    guides(size = "legend",color = "none") +
    ylim(0,0.004) -> plt
ggsave("../fig/WHR_LDL.pdf",plt,height = 6,width = 6.5)

WHRadjBMI_cs$label = ifelse(WHRadjBMI_cs$p_diff_HDL > -log10(0.05/nrow(WHRadjBMI_cs)), WHRadjBMI_cs$NearestGene, "")
ggplot(WHRadjBMI_cs, aes(x = diff_beta_HDL, y = diff_beta_WHRadjBMI, color = significance, size = p_diff_HDL, label = label)) +
    geom_point() +
    theme_bw() +
    theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
    xlab(expression(HDL:"|"~beta[F]-beta[M]~"|")) +
    ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
    labs(color = "", size = expression(-log[10](p[GxSex~HDL]))) +
    scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
    scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
    geom_text_repel(show.legend = F, max.overlaps = 50) +
    guides(size = "legend",color = "none") +
    ylim(0,0.004) -> plt
ggsave("../fig/WHR_HDL.pdf",plt,height = 6,width = 6.5)

WHRadjBMI_cs$label = ifelse(WHRadjBMI_cs$p_diff_TG > -log10(0.05/nrow(WHRadjBMI_cs)), WHRadjBMI_cs$NearestGene, "")
ggplot(WHRadjBMI_cs, aes(x = diff_beta_TG, y = diff_beta_WHRadjBMI, color = significance, size = p_diff_TG, label = label)) +
    geom_point() +
    theme_bw() +
    theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
    xlab(expression(TG:"|"~beta[F]-beta[M]~"|")) +
    ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
    labs(color = "", size = expression(-log[10](p[GxSex~TG]))) +
    scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
    scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
    geom_text_repel(show.legend = F, max.overlaps = 50) +
    guides(size = "legend",color = "none") +
    ylim(0,0.004) -> plt
ggsave("../fig/WHR_TG.pdf",plt,height = 6,width = 6.5)

WHRadjBMI_cs$label = ifelse(WHRadjBMI_cs$p_diff_maf > -log10(0.05/nrow(WHRadjBMI_cs)), WHRadjBMI_cs$NearestGene, "")
ggplot(WHRadjBMI_cs, aes(x = diff_maf, y = diff_beta_WHRadjBMI, color = significance, size = p_diff_maf, label = label)) +
    geom_point() +
    theme_bw() +
    theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
    xlab(expression("|"~MAF[F]-MAF[M]~"|")) +
    ylab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
    labs(color = "", size = expression(-log[10](p[sex~diff.~MAF]))) +
    scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
    scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
    geom_text_repel(show.legend = F, max.overlaps = 50) +
    guides(size = "legend",color = "none") +
    ylim(0,0.004) -> plt
ggsave("../fig/WHR_MAF.pdf",plt,height = 6,width = 6.5)

WHRadjBMI_cs$label = ifelse(WHRadjBMI_cs$p_diff_SHBG > -log10(0.05/nrow(WHRadjBMI_cs)), WHRadjBMI_cs$NearestGene, "")
ggplot(WHRadjBMI_cs, aes(x = diff_beta_WHRadjBMI, y = diff_beta_WHRadjBMI_adjshbg, color = significance, size = p_diff_SHBG, label = label)) +
    geom_abline(slope = 1, lty = 2, col = "darkgrey") +
    geom_point() +
    theme_bw() +
    theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
    xlab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
    ylab(expression(WHRadjBMI~adj.~SHBG:"|"~beta[F]-beta[M]~"|")) +
    labs(color = "", size = expression(-log[10](p[GxSex~SHBG]))) +
    scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
    scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
    geom_text_repel(show.legend = F, max.overlaps = 50) +
    guides(color = "none") -> plt
ggsave("../fig/WHRadjshbg.pdf",plt,height = 6,width = 6.5)

WHRadjBMI_cs$label = ifelse(WHRadjBMI_cs$p_diff_LDL > -log10(0.05/nrow(WHRadjBMI_cs)), WHRadjBMI_cs$NearestGene, "")
ggplot(WHRadjBMI_cs, aes(x = diff_beta_WHRadjBMI, y = diff_beta_WHRadjBMI_adjlipid, size = p_diff_LDL, color = significance, label = label)) +
    geom_abline(slope = 1, lty = 2, col = "darkgrey") +
    geom_point() +
    theme_bw() +
    theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
    xlab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
    ylab(expression(WHRadjBMI~adj.~lipids:"|"~beta[F]-beta[M]~"|")) +
    labs(color = "", size = expression(-log[10](p[GxSex~LDL]))) +
    scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
    scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
    geom_text_repel(show.legend = F, max.overlaps = 50) +
    guides(color = "none") -> plt
ggsave("../fig/WHRadjlipid_LDL.pdf",plt,height = 6,width = 6.5)

WHRadjBMI_cs$label = ifelse(WHRadjBMI_cs$p_diff_TG > -log10(0.05/nrow(WHRadjBMI_cs)), WHRadjBMI_cs$NearestGene, "")
ggplot(WHRadjBMI_cs, aes(x = diff_beta_WHRadjBMI, y = diff_beta_WHRadjBMI_adjlipid, size = p_diff_TG, color = significance, label = label)) +
    geom_abline(slope = 1, lty = 2, col = "darkgrey") +
    geom_point() +
    theme_bw() +
    theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
    xlab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
    ylab(expression(WHRadjBMI~adj.~lipids:"|"~beta[F]-beta[M]~"|")) +
    labs(color = "", size = expression(-log[10](p[GxSex~TG]))) +
    scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
    scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
    geom_text_repel(show.legend = F, max.overlaps = 50) +
    guides(color = "none") -> plt
ggsave("../fig/WHRadjlipid_TG.pdf",plt,height = 6,width = 6.5)

WHRadjBMI_cs$label = ifelse(WHRadjBMI_cs$p_diff_HDL > -log10(0.05/nrow(WHRadjBMI_cs)), WHRadjBMI_cs$NearestGene, "")
ggplot(WHRadjBMI_cs, aes(x = diff_beta_WHRadjBMI, y = diff_beta_WHRadjBMI_adjlipid, size = p_diff_HDL, color = significance, label = label)) +
    geom_abline(slope = 1, lty = 2, col = "darkgrey") +
    geom_point() +
    theme_bw() +
    theme(axis.title = element_text(size = 21),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
    xlab(expression(WHRadjBMI:"|"~beta[F]-beta[M]~"|")) +
    ylab(expression(WHRadjBMI~adj.~lipids:"|"~beta[F]-beta[M]~"|")) +
    labs(color = "", size = expression(-log[10](p[GxSex~HDL]))) +
    scale_color_manual(values = c("red","mediumvioletred","darkgrey")) +
    scale_size_continuous(labels = c(2,6,10),limits = c(0,10),breaks = c(2,6,10)) +
    geom_text_repel(show.legend = F, max.overlaps = 50) +
    guides(color = "none") -> plt
ggsave("../fig/WHRadjlipid_HDL.pdf",plt,height = 6,width = 6.5)