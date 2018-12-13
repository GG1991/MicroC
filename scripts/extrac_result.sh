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
		${AWK} '/n =/{print $3 "\t" $6"\t" $9}' $file >> tvsn_${pc}.txt
		tot_time=$(${AWK} '/n =/{print $6}' $file)
		tot_its=$(${AWK} '/n =/{print $9}' $file)

		${AWK} -v tot_time=${tot_time} -v tot_its=${tot_its} '
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
				printf("%-4d\t%e\t%e\t%e\t%e\n",
				$1, ($1 * tot_time / tot_its), $6, $10, $12);
			}

		}' $file > rvsi_${pc}_${n}.txt

	done
done
