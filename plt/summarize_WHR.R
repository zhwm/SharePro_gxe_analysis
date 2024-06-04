library(data.table)

get_res = function(rdir){
    data_lst = list()
    for (i in 1:22) {
        spath = paste0(rdir, i, '.sharepro.gxe.txt')
        data_lst[[i]] = fread(spath)
    }
    df = do.call(rbind, data_lst)
}

get_fa = function(rdir){
    data_lst = list()
    for (i in 1:22) {
        spath = paste0(rdir, i, '.fastGWA')
        data_lst[[i]] = fread(spath)
    }
    df = do.call(rbind, data_lst)
}

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

WHRadjBMI = get_res('../dat/WHRadjBMI_sex_FM/WHRadjBMI_FM_')
WHRadjBMI$nlog10pGxE_0_1 = (pchisq((WHRadjBMI$BETA0-WHRadjBMI$BETA1)**2/(WHRadjBMI$SE0^2+WHRadjBMI$SE1^2), df=1, lower.tail = F,log.p = T)/-log(10))
WHRadjBMI_cs = WHRadjBMI[WHRadjBMI$cs!=0]
WHRadjBMI_cs$NearestGene = mapply(get_nearest_gene, 
                        chr = WHRadjBMI_cs$CHR, 
                        pos = WHRadjBMI_cs$POS, 
                        MoreArgs = list(gtf = gtf))
WHRadjBMI_cs$idx1 = paste0(WHRadjBMI_cs$CHR, ".", WHRadjBMI_cs$POS, ".", WHRadjBMI_cs$A1, ".", WHRadjBMI_cs$A2)
WHRadjBMI_cs$idx2 = paste0(WHRadjBMI_cs$CHR, ".", WHRadjBMI_cs$POS, ".", WHRadjBMI_cs$A2, ".", WHRadjBMI_cs$A1)
WHRadjBMI_cs = as.data.frame(WHRadjBMI_cs)

ss = as.data.frame(fread('../dat/WHRadjBMI_sex_FM/GCST90012106_buildGRCh37.tsv.gz'))
ss$idx = paste0(ss$chromosome, ".", ss$base_pair_location, ".", ss$effect_allele, ".", ss$other_allele)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_SHBG_female = NA
WHRadjBMI_cs$se_SHBG_female = NA
WHRadjBMI_cs$af_SHBG_female = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_SHBG_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$beta
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_SHBG_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$standard_error
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_SHBG_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$effect_allele_frequency
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_SHBG_female = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$beta
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_SHBG_female = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$standard_error
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_SHBG_female = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$effect_allele_frequency

ss = as.data.frame(fread('../dat/WHRadjBMI_sex_FM/GCST90012108_buildGRCh37.tsv.gz'))
ss$idx = paste0(ss$chromosome, ".", ss$base_pair_location, ".", ss$effect_allele, ".", ss$other_allele)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_SHBG_male = NA
WHRadjBMI_cs$se_SHBG_male = NA
WHRadjBMI_cs$af_SHBG_male = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_SHBG_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$beta
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_SHBG_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$standard_error
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_SHBG_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$effect_allele_frequency
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_SHBG_male = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$beta
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_SHBG_male = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$standard_error
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_SHBG_male = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$effect_allele_frequency

ss = as.data.frame(fread('../dat/WHRadjBMI_sex_FM/GCST90012102_buildGRCh37.tsv.gz'))
ss$idx = paste0(ss$chromosome, ".", ss$base_pair_location, ".", ss$effect_allele, ".", ss$other_allele)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_T_female = NA
WHRadjBMI_cs$se_T_female = NA
WHRadjBMI_cs$af_T_female = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_T_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$beta
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_T_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$standard_error
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_T_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$effect_allele_frequency
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_T_female = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$beta
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_T_female = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$standard_error
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_T_female = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$effect_allele_frequency

ss = as.data.frame(fread('../dat/WHRadjBMI_sex_FM/GCST90012103_buildGRCh37.tsv.gz'))
ss$idx = paste0(ss$chromosome, ".", ss$base_pair_location, ".", ss$effect_allele, ".", ss$other_allele)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_T_male = NA
WHRadjBMI_cs$se_T_male = NA
WHRadjBMI_cs$af_T_male = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_T_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$beta
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_T_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$standard_error
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_T_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$effect_allele_frequency
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_T_male = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$beta
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_T_male = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$standard_error
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_T_male = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$effect_allele_frequency

ss = as.data.frame(fread('../dat/WHRadjBMI_sex_FM/HDL_FE_INT.fastGWA'))
ss$idx = paste0(ss$CHR, ".", ss$POS, ".", ss$A1, ".", ss$A2)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_HDL_female = NA
WHRadjBMI_cs$se_HDL_female = NA
WHRadjBMI_cs$af_HDL_female = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_HDL_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_HDL_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_HDL_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$AF1
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_HDL_female = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_HDL_female = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_HDL_female = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$AF1

ss = as.data.frame(fread('../dat/WHRadjBMI_sex_FM/HDL_MA_INT.fastGWA'))
ss$idx = paste0(ss$CHR, ".", ss$POS, ".", ss$A1, ".", ss$A2)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_HDL_male = NA
WHRadjBMI_cs$se_HDL_male = NA
WHRadjBMI_cs$af_HDL_male = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_HDL_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_HDL_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_HDL_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$AF1
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_HDL_male = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_HDL_male = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_HDL_male = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$AF1

ss = as.data.frame(fread('../dat/WHRadjBMI_sex_FM/LDL_FE_INT.fastGWA'))
ss$idx = paste0(ss$CHR, ".", ss$POS, ".", ss$A1, ".", ss$A2)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_LDL_female = NA
WHRadjBMI_cs$se_LDL_female = NA
WHRadjBMI_cs$af_LDL_female = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_LDL_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_LDL_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_LDL_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$AF1
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_LDL_female = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_LDL_female = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_LDL_female = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$AF1

ss = as.data.frame(fread('../dat/WHRadjBMI_sex_FM/LDL_MA_INT.fastGWA'))
ss$idx = paste0(ss$CHR, ".", ss$POS, ".", ss$A1, ".", ss$A2)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_LDL_male = NA
WHRadjBMI_cs$se_LDL_male = NA
WHRadjBMI_cs$af_LDL_male = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_LDL_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_LDL_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_LDL_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$AF1
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_LDL_male = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_LDL_male = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_LDL_male = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$AF1

ss = as.data.frame(fread('../dat/WHRadjBMI_sex_FM/TG_FE_INT.fastGWA'))
ss$idx = paste0(ss$CHR, ".", ss$POS, ".", ss$A1, ".", ss$A2)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_TG_female = NA
WHRadjBMI_cs$se_TG_female = NA
WHRadjBMI_cs$af_TG_female = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_TG_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_TG_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_TG_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$AF1
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_TG_female = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_TG_female = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_TG_female = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$AF1

ss = as.data.frame(fread('../dat/WHRadjBMI_sex_FM/TG_MA_INT.fastGWA'))
ss$idx = paste0(ss$CHR, ".", ss$POS, ".", ss$A1, ".", ss$A2)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_TG_male = NA
WHRadjBMI_cs$se_TG_male = NA
WHRadjBMI_cs$af_TG_male = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_TG_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_TG_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_TG_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$AF1
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_TG_male = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_TG_male = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_TG_male = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$AF1

ss = as.data.frame(get_fa('../dat/WHRadjBMI_sex_FM/WHRadjBMIshbg_FE_'))
ss$idx = paste0(ss$CHR, ".", ss$POS, ".", ss$A1, ".", ss$A2)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_adjshbg_female = NA
WHRadjBMI_cs$se_adjshbg_female = NA
WHRadjBMI_cs$af_adjshbg_female = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_adjshbg_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_adjshbg_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_adjshbg_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$AF1
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_adjshbg_female = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_adjshbg_female = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_adjshbg_female = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$AF1

ss = as.data.frame(get_fa('../dat/WHRadjBMI_sex_FM/WHRadjBMIshbg_MA_'))
ss$idx = paste0(ss$CHR, ".", ss$POS, ".", ss$A1, ".", ss$A2)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_adjshbg_male = NA
WHRadjBMI_cs$se_adjshbg_male = NA
WHRadjBMI_cs$af_adjshbg_male = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_adjshbg_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_adjshbg_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_adjshbg_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$AF1
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_adjshbg_male = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_adjshbg_male = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_adjshbg_male = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$AF1

ss = as.data.frame(get_fa('../dat/WHRadjBMI_sex_FM/WHRadjBMIlipid_FE_'))
ss$idx = paste0(ss$CHR, ".", ss$POS, ".", ss$A1, ".", ss$A2)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_adjlipid_female = NA
WHRadjBMI_cs$se_adjlipid_female = NA
WHRadjBMI_cs$af_adjlipid_female = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_adjlipid_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_adjlipid_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_adjlipid_female = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$AF1
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_adjlipid_female = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_adjlipid_female = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_adjlipid_female = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$AF1

ss = as.data.frame(get_fa('../dat/WHRadjBMI_sex_FM/WHRadjBMIlipid_MA_'))
ss$idx = paste0(ss$CHR, ".", ss$POS, ".", ss$A1, ".", ss$A2)
rownames(ss) = ss$idx
WHRadjBMI_cs$beta_adjlipid_male = NA
WHRadjBMI_cs$se_adjlipid_male = NA
WHRadjBMI_cs$af_adjlipid_male = NA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$beta_adjlipid_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$se_adjlipid_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx1 %in% rownames(ss),]$af_adjlipid_male = ss[WHRadjBMI_cs$idx1[WHRadjBMI_cs$idx1 %in% rownames(ss)],]$AF1
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$beta_adjlipid_male = -1 * ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$BETA
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$se_adjlipid_male = ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$SE
WHRadjBMI_cs[WHRadjBMI_cs$idx2 %in% rownames(ss),]$af_adjlipid_male = 1 - ss[WHRadjBMI_cs$idx2[WHRadjBMI_cs$idx2 %in% rownames(ss)],]$AF1

WHRadjBMI_cs$AF = (WHRadjBMI_cs$AF10 * WHRadjBMI_cs$N0 + WHRadjBMI_cs$AF11 * WHRadjBMI_cs$N1) / (WHRadjBMI_cs$N0 + WHRadjBMI_cs$N1)
WHRadjBMI_cs$significance = ifelse(WHRadjBMI_cs$nlog10pGxE_0_1 > -log10(5e-8), "Genome-wide GxSex in WHRadjBMI",ifelse(WHRadjBMI_cs$nlog10pGxE_0_1 > -log10(0.05 / nrow(WHRadjBMI_cs)), "Effect group-level GxSex in WHRadjBMI", "Non-significant"))
WHRadjBMI_cs$significance = factor(WHRadjBMI_cs$significance, levels = c("Genome-wide GxSex in WHRadjBMI","Effect group-level GxSex in WHRadjBMI","Non-significant"))
WHRadjBMI_cs$p_diff_SHBG = -log10(pchisq(((WHRadjBMI_cs$beta_SHBG_female - WHRadjBMI_cs$beta_SHBG_male) / sqrt(WHRadjBMI_cs$se_SHBG_female^2 + WHRadjBMI_cs$se_SHBG_male^2))^2, df = 1, lower.tail = F))
WHRadjBMI_cs$p_diff_T = -log10(pchisq(((WHRadjBMI_cs$beta_T_female - WHRadjBMI_cs$beta_T_male) / sqrt(WHRadjBMI_cs$se_T_female^2 + WHRadjBMI_cs$se_T_male^2))^2, df = 1, lower.tail = F))
WHRadjBMI_cs$p_diff_HDL = -log10(pchisq(((WHRadjBMI_cs$beta_HDL_female - WHRadjBMI_cs$beta_HDL_male) / sqrt(WHRadjBMI_cs$se_HDL_female^2 + WHRadjBMI_cs$se_HDL_male^2))^2, df = 1, lower.tail = F))
WHRadjBMI_cs$p_diff_LDL = -log10(pchisq(((WHRadjBMI_cs$beta_LDL_female - WHRadjBMI_cs$beta_LDL_male) / sqrt(WHRadjBMI_cs$se_LDL_female^2 + WHRadjBMI_cs$se_LDL_male^2))^2, df = 1, lower.tail = F))
WHRadjBMI_cs$p_diff_TG = -log10(pchisq(((WHRadjBMI_cs$beta_TG_female - WHRadjBMI_cs$beta_TG_male) / sqrt(WHRadjBMI_cs$se_TG_female^2 + WHRadjBMI_cs$se_TG_male^2))^2, df = 1, lower.tail = F))
WHRadjBMI_cs$p_diff_maf = -log10(pchisq(( (WHRadjBMI_cs$AF10 - WHRadjBMI_cs$AF11) / sqrt( WHRadjBMI_cs$AF * (1 - WHRadjBMI_cs$AF) * 0.5 * (1 / WHRadjBMI_cs$N0 + 1 / WHRadjBMI_cs$N1 ) ) )^2, df = 1, lower.tail = F))

WHRadjBMI_cs$p_diff_SHBG[WHRadjBMI_cs$p_diff_SHBG > 10] = 10
WHRadjBMI_cs$p_diff_T[WHRadjBMI_cs$p_diff_T > 10] = 10
WHRadjBMI_cs$p_diff_HDL[WHRadjBMI_cs$p_diff_HDL > 10] = 10
WHRadjBMI_cs$p_diff_LDL[WHRadjBMI_cs$p_diff_LDL > 10] = 10
WHRadjBMI_cs$p_diff_TG[WHRadjBMI_cs$p_diff_TG > 10] = 10

WHRadjBMI_cs$diff_beta_WHRadjBMI = abs(WHRadjBMI_cs$BETA0 * sqrt(2 * WHRadjBMI_cs$AF10 * (1 - WHRadjBMI_cs$AF10)) - WHRadjBMI_cs$BETA1 * sqrt(2 * WHRadjBMI_cs$AF11 * (1 - WHRadjBMI_cs$AF11)))
WHRadjBMI_cs$diff_beta_SHBG = abs(WHRadjBMI_cs$beta_SHBG_female * sqrt(2 * WHRadjBMI_cs$af_SHBG_female * (1 - WHRadjBMI_cs$af_SHBG_female)) - WHRadjBMI_cs$beta_SHBG_male * sqrt(2 * WHRadjBMI_cs$af_SHBG_male * (1 - WHRadjBMI_cs$af_SHBG_male)))
WHRadjBMI_cs$diff_beta_T = abs(WHRadjBMI_cs$beta_T_female * sqrt(2 * WHRadjBMI_cs$af_T_female * (1 - WHRadjBMI_cs$af_T_female)) - WHRadjBMI_cs$beta_T_male * sqrt(2 * WHRadjBMI_cs$af_T_male * (1 - WHRadjBMI_cs$af_T_male)))
WHRadjBMI_cs$diff_beta_HDL = abs(WHRadjBMI_cs$beta_HDL_female * sqrt(2 * WHRadjBMI_cs$af_HDL_female * (1 - WHRadjBMI_cs$af_HDL_female)) - WHRadjBMI_cs$beta_HDL_male * sqrt(2 * WHRadjBMI_cs$af_HDL_male * (1 - WHRadjBMI_cs$af_HDL_male)))
WHRadjBMI_cs$diff_beta_LDL = abs(WHRadjBMI_cs$beta_LDL_female * sqrt(2 * WHRadjBMI_cs$af_LDL_female * (1 - WHRadjBMI_cs$af_LDL_female)) - WHRadjBMI_cs$beta_LDL_male * sqrt(2 * WHRadjBMI_cs$af_LDL_male * (1 - WHRadjBMI_cs$af_LDL_male)))
WHRadjBMI_cs$diff_beta_TG = abs(WHRadjBMI_cs$beta_TG_female * sqrt(2 * WHRadjBMI_cs$af_TG_female * (1 - WHRadjBMI_cs$af_TG_female)) - WHRadjBMI_cs$beta_TG_male * sqrt(2 * WHRadjBMI_cs$af_TG_male * (1 - WHRadjBMI_cs$af_TG_male)))
WHRadjBMI_cs$diff_maf = abs(WHRadjBMI_cs$AF10 - WHRadjBMI_cs$AF11)
WHRadjBMI_cs$diff_beta_WHRadjBMI_adjshbg = abs(WHRadjBMI_cs$beta_adjshbg_female * sqrt(2 * WHRadjBMI_cs$AF10 * (1 - WHRadjBMI_cs$AF10)) - WHRadjBMI_cs$beta_adjshbg_male * sqrt(2 * WHRadjBMI_cs$AF11 * (1 - WHRadjBMI_cs$AF11)))
WHRadjBMI_cs$diff_beta_WHRadjBMI_adjlipid = abs(WHRadjBMI_cs$beta_adjlipid_female * sqrt(2 * WHRadjBMI_cs$AF10 * (1 - WHRadjBMI_cs$AF10)) - WHRadjBMI_cs$beta_adjlipid_male * sqrt(2 * WHRadjBMI_cs$AF11 * (1 - WHRadjBMI_cs$AF11)))

write.table(WHRadjBMI_cs, file = '../doc/WHR_summary.txt', quote=F, sep='\t', row.names=F, col.names=T)