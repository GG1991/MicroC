#!/bin/bash

source ./vars.sh

mkdir -p result

export OMP_NUM_THREADS=1

PATH="./result"
EXEC="../build/test/test_cg"

for pc in ${petsc_pc[@]}; do
	for n in ${sizes[@]}; do
		for a in ${factor_a[@]}; do

			echo "${EXEC} $n $a -pc_type $pc"
			${EXEC} $n $a -pc_type $pc -ksp_monitor > "$PATH/out_${pc}_${n}_${a}.txt"


		done
	done
done
