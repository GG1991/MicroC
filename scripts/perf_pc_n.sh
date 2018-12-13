#!/bin/bash

petsc_pc=( "jacobi" "ilu" "mg" "sor" "gasm" ) 
sizes=( 10 15 20 25 30 ) 

mkdir -p result

PATH="./result"
EXEC="../build/test/test_solver_linear"

for pc in ${petsc_pc[@]}; do
	for n in ${sizes[@]}; do

		echo "${EXEC} $n -pc_type $pc"
		${EXEC} $n -pc_type $pc -ksp_monitor > "$PATH/out_${pc}_${n}.txt"

	done
done