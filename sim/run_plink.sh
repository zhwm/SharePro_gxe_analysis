module load plink/1.9b_6.21-x86_64
cd ${1}
cd ${2}

awk '{print $1"\t"$2"\t"$3+1}' L.phen > L.cond
awk '{print $1"\t"$2"\t"$3+1}' M.phen > M.cond
awk '{print $1"\t"$2"\t"$3+1}' S.phen > S.cond

st2=`date +%s.%N`
for i in $(seq 1 50); do plink --bfile ../${1} --gxe --covar L.cond --out CL$i --pheno L.phen --pheno-name trait$i; done
for i in $(seq 1 50); do plink --bfile ../${1} --gxe --covar M.cond --out CM$i --pheno M.phen --pheno-name trait$i; done
for i in $(seq 1 50); do plink --bfile ../${1} --gxe --covar S.cond --out CS$i --pheno S.phen --pheno-name trait$i; done
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > plink.time
