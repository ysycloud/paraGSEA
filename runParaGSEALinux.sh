#!/bin/bash

#param list :filename topn
#quick_search_serial -i data/data_for_test.txt -n 10

#param list :filename thread_num topn
#quick_search_omp -i data/data_for_test.txt -t 4 -n 10 

#param list :process_num pernum hostfile filename topn
#mpirun -n 2 -ppn 2 -hostfile example/hostfile quick_search_mpi -i data/data_for_test.txt -n 15

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile ES_Matrix_ompi_nocom -t 4 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix_test

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile ES_Matrix_ompi_p2p -t 4 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix_test

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile ES_Matrix_ompi_cocom -t 4 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix_test

#param list :process_num pernum hostfile thread_num cluster_num filename outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile Cluster_KMediods_ompi -t 4 -c 5 -i data/ES_Matrix_test -o data/Cluster_result_test.txt

#param list :process_num pernum hostfile thread_num cluster_num filename outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile Cluster_KMediods++_ompi -t 4 -c 5 -i data/ES_Matrix_test -o data/Cluster_result_test.txt