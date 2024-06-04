module load rstudio-server
cd ${1}
mkdir ${2}\_${3}\_${4}
cd ${2}\_${3}\_${4}
Rscript ../../simulate_traits.R ../${1} ${2} 0 ${3} ${4} 0
