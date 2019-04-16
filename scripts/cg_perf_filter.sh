#!/bin/bash

source vars.sh

mkdir -p result

PATH="./result"
AWK="/usr/bin/awk"

for ksp in ${ksp_type[@]}; do
for pc in ${pc_type[@]}; do
for n in ${sizes[@]}; do
for tol in ${tolerances[@]}; do
for a in ${factor_a[@]}; do

	filein="$PATH/out_${ksp}_${pc}_${n}_${a}_${tol}.txt"
	fileout="$PATH/out_${ksp}_${pc}_${n}_${a}_${tol}_filtered.txt"
	echo "Its      Time     RNorm    Norm" > $fileout

	solver_time=$(${AWK} '/solve_v/{print $3}' $filein)
	tot_its=$(${AWK} '/n =/{print $9}' $filein)

	${AWK} -v solver_time=${solver_time} -v tot_its=${tot_its} '
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
			printf("%-4d\t%e\t%e\t%e\n",
			$1, ($1 * solver_time / tot_its), $6, $10);
		}

	}' $filein >> $fileout

done
done
done
done
done
