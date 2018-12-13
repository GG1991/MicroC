#!/bin/bash

petsc_pc=( "jacobi" "ilu" "mg" "sor" "gasm" )
sizes=( 10 15 20 25 30 )

mkdir -p result

PATH="./result"

for pc in ${petsc_pc[@]}; do
	echo "#   N     time" > tvsn_${pc}.txt
	for n in ${sizes[@]}; do

		echo "Extract from $PATH/out_${pc}_${n}.txt"
		/usr/bin/awk '/n =/{print $3 "\t" $6}' "$PATH/out_${pc}_${n}.txt" >> tvsn_${pc}.txt

	done
done
