#!/bin/bash

#param list :filename topn
#./bin/quick_search_serial "data/data_for_test.txt" 10

#param list :filename thread_num topn paraway(0/1)
#./bin/quick_search_omp "data/data_for_test.txt" 4 10 1

#param list :process_num pernum hostfile filename topn
#mpirun -n 2 -ppn 2 -hostfile hostfile ./bin/quick_search_mpi "data/data_for_test.txt" 15

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
#mpirun -n 2 -ppn 2 -hostfile hostfile ./bin/ES_Matrix_ompi_nocom 4 50 "data/data_for_test.txt" "data/data_for_test.txt" "data/ES_Matrix_test"

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
#mpirun -n 2 -ppn 2 -hostfile hostfile ./bin/ES_Matrix_ompi_p2p 4 50 "data/data_for_test.txt" "data/data_for_test.txt" "data/ES_Matrix_test"

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
#mpirun -n 2 -ppn 2 -hostfile hostfile ./bin/ES_Matrix_ompi_cocom 4 50 "data/data_for_test.txt" "data/data_for_test.txt" "data/ES_Matrix_test"

#param list :process_num pernum hostfile thread_num cluster_num filename outfilename
#mpirun -n 2 -ppn 2 -hostfile hostfile ./bin/Cluster_KMediods_ompi 4 12 "data/ES_Matrix_test" "data/Cluster_result_test.txt"

#param list :process_num pernum hostfile thread_num cluster_num filename outfilename
#mpirun -n 2 -ppn 2 -hostfile hostfile ./bin/Cluster_KMediods++_ompi 4 12 "data/ES_Matrix_test" "data/Cluster_result_test.txt"