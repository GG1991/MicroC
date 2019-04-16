#!/bin/bash

source ./vars.sh

mkdir -p result

export OMP_NUM_THREADS=1

OPATH="./result"
EXEC="../build/test/test_cg"

for ksp in ${ksp_type[@]}; do
for pc in ${pc_type[@]}; do
for n in ${sizes[@]}; do
for tol in ${tolerances[@]}; do
for a in ${factor_a[@]}; do

	echo "${EXEC} $n $a -ksp_type ${ksp} -pc_type ${pc} -ksp_rtol $tol -ksp_monitor -kps_rtol $tol -ksp_monitor_true_residual" \
		| tee "$OPATH/out_${ksp}_${pc}_${n}_${a}_${tol}.txt"
	${EXEC} $n $a -ksp_type ${ksp} -pc_type ${pc} -ksp_rtol $tol \
		-ksp_monitor_true_residual \
		>> "$OPATH/out_${ksp}_${pc}_${n}_${a}_${tol}.txt"

done
done
done
done
done
