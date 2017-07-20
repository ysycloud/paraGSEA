#!/bin/bash
read -p "Please enter the way to execute the quick_search(0_serial,1_openmp,2_mpi,3_profile):" quick_search_way

# execute matlab script to parse the data(pre_treatment)
runPreGSEAbyMatlab.sh
#runparaPreGSEAbyMatlab.sh

#quick_search 
case "$quick_search_way" in
	0)quick_search_serial -i ../data/data_for_test.txt -n 8 -s ../data/data_for_test_cidnum.txt -r ../data/Reference;;
	1)quick_search_omp -i ../data/data_for_test.txt -t 4 -n 10 -s ../data/data_for_test_cidnum.txt -r ../data/Reference;;
	2)mpirun -n 2 -ppn 2 -hostfile hostfile quick_search_mpi -i ../data/data_for_test.txt -n 8 -s ../data/data_for_test_cidnum.txt -r ../data/Reference
	3)mpirun -n 2 -ppn 2 -hostfile hostfile quick_search_profile -i ../data/data_for_test.txt -n 8 -t 4 -l 50 -s ../data/data_for_test_cidnum.txt -r ../data/Reference
esac