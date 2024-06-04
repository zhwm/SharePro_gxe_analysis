module load StdEnv/2020 r/4.0.0
cd ${1}
cd ${2}

for i in $(seq 1 50); do Rscript ../../susie.R CL$i ../${1}.ld; done
for i in $(seq 1 50); do Rscript ../../susie.R CM$i ../${1}.ld; done
for i in $(seq 1 50); do Rscript ../../susie.R CS$i ../${1}.ld; done


