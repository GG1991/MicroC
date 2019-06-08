#!/bin/bash

#pc_type=( "jacobi" ) 
#ksp_type=( "cg" "bcgs" "gmres" "bicg" ) 
pc_type=( "none" "jacobi" "ilu" "mg") 
ksp_type=( "cg" ) 
sizes=( 50 )

mkdir -p result

export OMP_NUM_THREADS=1

OPATH="./result"

for ksp in ${ksp_type[@]}; do
for pc in ${pc_type[@]}; do
for n in ${sizes[@]}; do

EXEC="../build/test/test_cg $n 10.0 -ksp_type ${ksp} -pc_type ${pc} -ksp_monitor_true_residual -ksp_rtol 1.0e-5 -ksp_atol 1.0e-20"

	fileo="$OPATH/out_${ksp}_${pc}_${n}.txt"

	echo "${EXEC}" | tee $fileo

	${EXEC}  > $fileo

	#file2="times-$ksp-$pc-homogeneous.txt"
	file2="times-$ksp-$pc-heterogeneous.txt"
	tim=$(awk '/time/{print $6}' $fileo)
	its=$(awk '/time/{print $9}' $fileo)
	awk -v its=$its -v tim=$tim '/precon/{print $1 "\t" (tim * $1 / its) "\t" $6 "\t" $10 "\t" $12}' $fileo > $file2
done
done
done
