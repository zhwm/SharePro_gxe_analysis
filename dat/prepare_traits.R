INTfun <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
setwd("~/scratch/SharePro_gxe/dat/GxE/")

library(data.table)

### load EUR data
traits = fread("trait.ukb.txt")
EUR = fread('~/scratch/UKB_Geno/all_idpEUR.fam')$V1
traits = traits[traits$f.eid %in% EUR,]

### process phenotypes
traits$Age <- traits$f.21022.0.0
traits$Age2 <- traits$f.21022.0.0^2
traits$centre <- as.factor(traits$f.54.0.0)
traits$array <- ifelse(traits$f.22000.0.0>0,1,0)
traits$FFR <- traits$f.20258.0.0
traits$BMI <- traits$f.21001.0.0
traits$WHR <- traits$f.48.0.0/traits$f.49.0.0
traits$smoke <- traits$f.20116.0.0
traits$sex <- traits$f.22001.0.0
traits$LDL <- traits$f.30780.0.0
traits$TG <- traits$f.30870.0.0
traits$HDL <- traits$f.30760.0.0

#bioavail_t <- function(ALB,SHBG,Testo) {
#  kat = 3.6e4
#  kt = 1e9
#  wt = 15046.6446
#  ALB = ALB * wt
#  a = kat + kt + (kat * kt) * (SHBG + ALB - Testo)
#  b = 1 + kt * SHBG + kat * ALB - (kat + kt) * Testo
#  free_Testo = (-b + sqrt(b^2 + 4 * a * Testo)) / 2 / a * 1e9
#  c = kat * ALB / wt / 69000
#  bioavail_Testo = free_Testo * (1 + c)
#  return(bioavail_Testo)
#}

### get INT smoking
get_INT <- function(df,trait,pop){
  df$resid <- NA
  df$resid[!is.na(df[[trait]])] <- 
    INTfun(summary(lm(df[[trait]] ~ Age + Age2 + array + centre + sex +
                        f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + 
                        f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10 + 
                        f.22009.0.11 + f.22009.0.12 + f.22009.0.13 + f.22009.0.14 + f.22009.0.15 + 
                        f.22009.0.16 + f.22009.0.17 + f.22009.0.18 + f.22009.0.19 + f.22009.0.20, 
                      data = df))$residuals)
  write.table(df[,c("f.eid","f.eid","resid")],
              paste0(trait,'_',pop,"_INT.txt"),sep=" ",col.names=F,row.names=F,quote=F)
}

### FFR vs Smoke

NS = traits[traits$smoke==0]
ES = traits[traits$smoke==1]
CS = traits[traits$smoke==2]

get_INT(NS,'FFR','NS')
get_INT(ES,'FFR','ES')
get_INT(CS,'FFR','CS')


### WHRadjBMI vs Sex

get_INT_sex_WHR_BMI <- function(df,trait,pop){
  df$resid <- NA
  df$resid[!is.na(df[[trait]]) & (!is.na(df$BMI))] <- 
    INTfun(summary(lm(df[[trait]] ~ Age + Age2 + array + centre + BMI +
                        f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + 
                        f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10 + 
                        f.22009.0.11 + f.22009.0.12 + f.22009.0.13 + f.22009.0.14 + f.22009.0.15 + 
                        f.22009.0.16 + f.22009.0.17 + f.22009.0.18 + f.22009.0.19 + f.22009.0.20, 
                      data = df))$residuals)
  write.table(df[,c("f.eid","f.eid","resid")],
              paste0(trait,'_',pop,"_INT.txt"),sep=" ",col.names=F,row.names=F,quote=F)
}

FE = traits[traits$sex==0]
MA = traits[traits$sex==1]

get_INT_sex_WHR_BMI(FE,'WHR','FE_BMI')
get_INT_sex_WHR_BMI(MA,'WHR','MA_BMI')


### biot vs Sex

get_INT_biot_sex <- function(df,trait,pop){
  df$resid <- NA
  df$resid[!is.na(df[[trait]])] <- 
    INTfun(summary(lm(df[[trait]] ~ Age + Age2 + array + centre +
                        f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + 
                        f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10 + 
                        f.22009.0.11 + f.22009.0.12 + f.22009.0.13 + f.22009.0.14 + f.22009.0.15 + 
                        f.22009.0.16 + f.22009.0.17 + f.22009.0.18 + f.22009.0.19 + f.22009.0.20, 
                      data = df))$residuals)
  write.table(df[,c("f.eid","f.eid","resid")],
              paste0(trait,'_',pop,"_INT.txt"),sep=" ",col.names=F,row.names=F,quote=F)
}

FE = traits[traits$sex==0]
MA = traits[traits$sex==1]

get_INT_biot_sex(FE,'TG','FE')
get_INT_biot_sex(MA,'TG','MA')

get_INT_biot_sex(FE,'LDL','FE')
get_INT_biot_sex(MA,'LDL','MA')

get_INT_biot_sex(FE,'HDL','FE')
get_INT_biot_sex(MA,'HDL','MA')

