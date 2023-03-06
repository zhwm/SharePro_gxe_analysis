module load plink/1.9b_6.21-x86_64
mkdir LOCI
cd LOCI
plink --bfile ~/scratch/UKB_Geno/CHR --chr CHR --from-bp ST --to-bp ED --keep ../UKB_A.fam --maf 0.01 --geno 0.01 --make-bed --out LOCI
plink --bfile LOCI --indep-pairwise 1000 80 0.1 --out LOCI
plink --bfile LOCI --extract LOCI.prune.in --recode A --out LOCI
st2=`date +%s.%N`
plink --bfile LOCI --matrix --r --out LOCI
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > getld.time
