source ~/py3/bin/activate
cd ${1}
cd ${2}

python ../../pre_SH.py --dir . --loci ${1}

st2=`date +%s.%N`
python ~/scratch/SharePro_gxe/src/sharepro_gxe.py --zld ACL.zld --zdir . --N 50000 --save SA --prefix ACL --verbose --K 10
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > ACL.time

st2=`date +%s.%N`
python ~/scratch/SharePro_gxe/src/sharepro_gxe.py --zld ACM.zld --zdir . --N 35000 --save SA --prefix ACM --verbose --K 10
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > ACM.time

st2=`date +%s.%N`
python ~/scratch/SharePro_gxe/src/sharepro_gxe.py --zld ACS.zld --zdir . --N 30000 --save SA --prefix ACS --verbose --K 10
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > ACS.time

st2=`date +%s.%N`
python ~/scratch/SharePro_gxe/src/sharepro_gxe.py --zld CS.zld --zdir . --N 25000 5000 --save SH --prefix CS --verbose --K 10 --sigma 1e-2
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > CS.time

st2=`date +%s.%N`
python ~/scratch/SharePro_gxe/src/sharepro_gxe.py --zld CM.zld --zdir . --N 25000 10000 --save SH --prefix CM --verbose --K 10 --sigma 1e-2
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > CM.time

st2=`date +%s.%N`
python ~/scratch/SharePro_gxe/src/sharepro_gxe.py --zld CL.zld --zdir . --N 25000 25000 --save SH --prefix CL --verbose --K 10 --sigma 1e-2
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > CL.time

rm *.z
rm *.log
