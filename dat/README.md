## GWAS

First, we performed self-reported smoking status stratified GWAS for FFR and sex stratified GWAS for WHRadjBMI in individuals of European ancestry in the UK Biobank.

```commandline
Rscript prepare_traits.R
~/utils/gcta_1.93.2beta/gcta64 --mbfile ~/scratch/UKB_Geno/geno_chrs_bgen_files.txt --fastGWA-lr --pheno SSSS --threads 40 --out SSSS --maf 0.001
```

## GxSmoking in FFR

We applied SharePro to these smoking status stratified GWAS summary statistic with the UK Biobank European ancestry LD matrix calculated by Weissbrod et al.

```commandline
python ~/scratch/SharePro_gxe/src/sharepro_gxe_ukb.py --ukb ~/utils/ukb/lst/GENO.lst --zdir FFR_CS_GENO.z FFR_ES_GENO.z FFR_NS_GENO.z --LDdir ~/scratch/UKBBLD/ ~/scratch/UKBBLD/ ~/scratch/UKBBLD/ --N $(cat FFR_CS.N) $(cat FFR_ES.N) $(cat FFR_NS.N) --save CEN --prefix FFR_smoke_CEN_GENO --K 5 --verbose
```

Additionally, we performed heterogeneity test with smoking status stratified GWAS summary statistics and obtained variant-level GxSmoking p-values.

```commandline
python FFR_smoke_gxe.py
```

We annotate each effect group to the nearest protein coding gene based on UCSC genome annotations and visualize results for GxSmoking analyses,

```commandline
python map2gene.py --gtf ~/utils/gencode.v19.protein_coding.gtf --rs ~/utils/ukb/idx/allrsid.dict --cs ../dat/GxE/FFR_smoke_CEN/CEN/FFR_smoke_CEN_ --save ../doc/FFR_smoke_CEN_cgene.txt --eQTL ../dat/anno/GTEx.eQTL.rsid.dict --sQTL ../dat/anno/GTEx.sQTL.rsid.dict
Rscript plot_FFR.R
```

## GxSex in fat distribution

Similarly, we investigated sex-differentiated genetic effects on body fat distribution, measured by waist-to-hip ratio adjusted for body mass index (WHRadjBMI).

```commandline
python ~/scratch/SharePro_gxe/src/sharepro_gxe_ukb.py --ukb ~/utils/ukb/lst/GENO.lst --zdir WHR_FE_BMI_GENO.z WHR_MA_BMI_GENO.z --LDdir ~/scratch/UKBBLD/ ~/scratch/UKBBLD/ --N $(cat WHR_FE_BMI.N) $(cat WHR_MA_BMI.N) --save FM --prefix WHR_BMI_FM_GENO --K 5 --verbose
```

Next, we investigated whether GxSex in WHRadjBMI may be explained by GxSex in sex hormones and lipid metabolism. 

```commandline
python map2gene.py --gtf ~/utils/gencode.v19.protein_coding.gtf --rs ~/utils/ukb/idx/allrsid.dict --cs ../dat/GxE/WHR_BMI_FM/FM/WHR_BMI_FM_ --save ../doc/WHR_cgene.txt --eQTL ../dat/anno/GTEx.eQTL.rsid.dict --sQTL ../dat/anno/GTEx.sQTL.rsid.dict
python matchss_WHR.py
```

Finally, we visualize the GxSex in fat distribution.
```commandline
Rscript plot_WHR.R
```