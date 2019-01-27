#!/bin/bash

source ./vars.sh

mkdir -p result

PATH="./result"
EXEC="../build/test/test_solver_linear"

for pc in ${petsc_pc[@]}; do
	for n in ${sizes[@]}; do
		for a in ${factor_a[@]}; do

			echo "${EXEC} $n -pc_type $pc"
			${EXEC} $n $a -pc_type $pc -ksp_monitor_true_residual > "$PATH/out_${pc}_${n}_${a}.txt"


		done
	done
done
