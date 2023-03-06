source ~/py3/bin/activate
cd ${1}
cd ${2}

st2=`date +%s.%N`
for i in $(seq 1 50); do ~/utils/GEM_1.4.4_static --pheno-file L.phen --delim \t --bfile ../${1} --sampleid-name IID --pheno-name trait$i --center 1 --exposure-name cond --out CL$i; done
for i in $(seq 1 50); do ~/utils/GEM_1.4.4_static --pheno-file M.phen --delim \t --bfile ../${1} --sampleid-name IID --pheno-name trait$i --center 1 --exposure-name cond --out CM$i; done
for i in $(seq 1 50); do ~/utils/GEM_1.4.4_static --pheno-file S.phen --delim \t --bfile ../${1} --sampleid-name IID --pheno-name trait$i --center 1 --exposure-name cond --out CS$i; done
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > gem.time
