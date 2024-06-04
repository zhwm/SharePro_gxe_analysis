module load StdEnv/2020 plink/1.9b_6.21-x86_64
cd ${1}
cd ${2}
for i in $(seq 1 50); do
  awk '{print$1"\t"$4"\t"$5"\t"$7"\t"$12"\t"$13"\t"$16"\t"$6}' CL$i | sed 1d | sed "1iSNP\tA1\tA2\tfreq\tb\tse\tP\tn" > CL$i\.ma;
  ~/utils/gcta_1.93.2beta/gcta64 --bfile ../${1} --cojo-file CL$i\.ma --cojo-slct --diff-freq 1.0 --out CL$i;
  plink --bfile ../${1} --clump CL$i\.ma --out CL$i --clump-p1 5e-8 --clump-kb 3000 --clump-r2 0.01;
done

for i in $(seq 1 50); do
  awk '{print$1"\t"$4"\t"$5"\t"$7"\t"$12"\t"$13"\t"$16"\t"$6}' CM$i | sed 1d | sed "1iSNP\tA1\tA2\tfreq\tb\tse\tP\tn" > CM$i\.ma;
  ~/utils/gcta_1.93.2beta/gcta64 --bfile ../${1} --cojo-file CM$i\.ma --cojo-slct --diff-freq 1.0 --out CM$i;
  plink --bfile ../${1} --clump CM$i\.ma --out CM$i --clump-p1 5e-8 --clump-kb 3000 --clump-r2 0.01;
done

for i in $(seq 1 50); do
  awk '{print$1"\t"$4"\t"$5"\t"$7"\t"$12"\t"$13"\t"$16"\t"$6}' CS$i | sed 1d | sed "1iSNP\tA1\tA2\tfreq\tb\tse\tP\tn" > CS$i\.ma;
  ~/utils/gcta_1.93.2beta/gcta64 --bfile ../${1} --cojo-file CS$i\.ma --cojo-slct --diff-freq 1.0 --out CS$i;
  plink --bfile ../${1} --clump CS$i\.ma --out CS$i --clump-p1 5e-8 --clump-kb 3000 --clump-r2 0.01;
done
