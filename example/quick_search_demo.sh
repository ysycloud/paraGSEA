#!/bin/bash
read -p "Please enter the way to execute the quick_search(0_serial,1_openmp,2_mpi):" quick_search_way

# execute matlab script to parse the data  
matlab -nodesktop -nosplash -nojvm -r "file_input='../data/modzs_n272x978.gctx'; file_name='../data/data_for_test.txt'; file_name_cid='../data/data_for_test_cid.txt'; PreGSEA; quit;"

#quick_search
case "$quick_search_way" in
	0)quick_search_serial "../data/data_for_test.txt" 10;;
	1)quick_search_omp "../data/data_for_test.txt" 4 10 1;;
	2)mpirun -n 2 -ppn 2 -hostfile hostfile quick_search_mpi "../data/data_for_test.txt" 10;;
esac