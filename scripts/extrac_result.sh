#!/bin/bash

source vars.sh

mkdir -p result

PATH="./result"
AWK="/usr/bin/awk"

for pc in ${petsc_pc[@]}; do
	echo "#   N     time" > tvsn_${pc}.txt
	for n in ${sizes[@]}; do

		file="$PATH/out_${pc}_${n}.txt"
		echo "Extract from $file"
		${AWK} '/n =/{print $3 "\t" $6}' "$file" >> tvsn_${pc}.txt

		${AWK} '
		BEGIN{flag = 0;}
		{
			if($1 == "SOLVER_START") {
				flag = 1;
				getline;
			}
			if($1=="SOLVER_END") {
				exit;
			}
			if(flag == 1) {
				print $1 "\t" $6 "\t" $10 "\t" $12;
			}

		}' $file > rvsi_${pc}_${n}.txt

	done
done
