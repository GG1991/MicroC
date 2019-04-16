#!/bin/bash

#pc_types=( "jacobi" "ilu" "mg") 
pc_type=( "jacobi" ) 
#ksp_type=( "cg" "bcgs" ) 
ksp_type=( "cg" ) 
#factor_a=( "1.0" "1000.0" )
factor_a=( "1.0" "100.0" "1000.0" )
tolerances=( "1.0e-10" )
#sizes=( 10 15 20 25 30)
sizes=( 70 )
