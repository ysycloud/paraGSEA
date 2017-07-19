#!/bin/bash

#param list : -i filename; -n topn; -s sample number file; -r reference directory  
#bin/quick_search_serial -i data/data_for_test.txt -n 8 -s data/data_for_test_cidnum.txt -r data/Reference

#param list : -i filename; -t thread_num; -n topn; -s sample number file; -r reference directory
#bin/quick_search_omp -i data/data_for_test.txt -t 4 -n 10 -s data/data_for_test_cidnum.txt -r data/Reference

#param list : -i filename; -n topn; -s sample number file; -r reference directory
#mpirun -n 2 -ppn 2 -hostfile example/hostfile bin/quick_search_mpi -i data/data_for_test.txt -n 8 -s data/data_for_test_cidnum.txt -r data/Reference

#param list : -i filename; -n topn; -t thread_num; -l siglen; -s sample number file; -r reference directory
#mpirun -n 2 -ppn 2 -hostfile example/hostfile bin/quick_search_profile -i data/data_for_test.txt -l 50 -t 4 -n 8 -s data/data_for_test_cidnum.txt -r data/Reference

#param list : -n process_num; -t thread_num; -l siglen; -1 filename1; -2 filename2; -p proportion; -w ifwrite; -o outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile bin/ES_Matrix_ompi_nocom -t 4 -l 50 -a 2 -p 1 -w 1 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix_test

#param list : -n process_num; -t thread_num; -l siglen; -1 filename1; -2 filename2; -p proportion; -w ifwrite; -o outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile bin/ES_Matrix_ompi_p2p -t 4 -l 50 -a 2 -p 1 -w 1 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix_test

#param list : -n process_num; -t thread_num; -l siglen; -1 filename1; -2 filename2; -p proportion; -w ifwrite; -o outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile bin/ES_Matrix_ompi_cocom -t 4 -l 50 -a 2 -p 1 -w 1 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix_test

#param list : -n process_num; -t thread_num; -c cluster_num; -w ifwrite; -i filename; -o outfilename; -s sample number file; -r reference directory 
#mpirun -n 2 -ppn 2 -hostfile example/hostfile bin/Cluster_KMediods_ompi -t 4 -c 5 -w 1 -i data/ES_Matrix_test -o data/Cluster_result_test.txt  -s data/data_for_test_cidnum.txt -r data/Reference

#param list : -n process_num; -t thread_num; -c cluster_num; -w ifwrite; -i filename; -o outfilename; -s sample number file; -r reference directory 
#mpirun -n 2 -ppn 2 -hostfile example/hostfile bin/Cluster_KMediods++_ompi -t 4 -c 5 -w 1 -i data/ES_Matrix_test -o data/Cluster_result_test.txt  -s data/data_for_test_cidnum.txt -r data/Reference
