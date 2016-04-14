#!/bin/bash

read -p "Please enter total number of processes:" -t 30 n
ppn=3
thread_num=4
siglen=50
cluster_num=12
inputfile="data/data_for_test.txt"
outputfile="data/Cluster_result_test.txt"

mpirun -n $n -ppn $ppn -hostfile hostfile ES_Matrix_ompi_cocom $thread_num $siglen $inputfile $inputfile "data/ES_Matrix_tmp"

mpirun -n $n -ppn $ppn -hostfile hostfile Cluster_KMediods++_ompi $thread_num $cluster_num "data/ES_Matrix_tmp" $outputfile


#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
#mpirun -n 2 -ppn 2 -hostfile hostfile ES_Matrix_ompi_nocom 4 50 "data/data_for_test.txt" "data/data_for_test.txt" "data/ES_Matrix_test"

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
#mpirun -n 2 -ppn 2 -hostfile hostfile ES_Matrix_ompi_p2p 4 50 "data/data_for_test.txt" "data/data_for_test.txt" "data/ES_Matrix_test"

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
#mpirun -n 2 -ppn 2 -hostfile hostfile ES_Matrix_ompi_cocom 4 50 "data/data_for_test.txt" "data/data_for_test.txt" "data/ES_Matrix_test"

#param list :process_num pernum hostfile thread_num cluster_num filename outfilename
#mpirun -n 2 -ppn 2 -hostfile hostfile Cluster_KMediods_ompi 4 12 "data/ES_Matrix_test" "data/Cluster_result_test.txt"

#param list :process_num pernum hostfile thread_num cluster_num filename outfilename
#mpirun -n 2 -ppn 2 -hostfile hostfile Cluster_KMediods++_ompi 4 12 "data/ES_Matrix_test" "data/Cluster_result_test.txt"