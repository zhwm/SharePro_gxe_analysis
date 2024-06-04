source ~/py3/bin/activate

cd ${1}
cd ${2}

for i in $(seq 1 50); do python ../../sharepro_gxe_sim.py --z CL$i\.qassoc.gxe --ld ../${1}.ld --save CL$i; done
for i in $(seq 1 50); do python ../../sharepro_gxe_sim.py --z CM$i\.qassoc.gxe --ld ../${1}.ld --save CM$i; done
for i in $(seq 1 50); do python ../../sharepro_gxe_sim.py --z CS$i\.qassoc.gxe --ld ../${1}.ld --save CS$i; done

for i in $(seq 1 50); do python ../../sharepro_gxe_combine.py --z CL$i --ld ../${1}.ld --save CL$i; done
for i in $(seq 1 50); do python ../../sharepro_gxe_combine.py --z CM$i --ld ../${1}.ld --save CM$i; done
for i in $(seq 1 50); do python ../../sharepro_gxe_combine.py --z CS$i --ld ../${1}.ld --save CS$i; done
