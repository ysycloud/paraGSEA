#!/bin/bash
read -p "Please enter the way to calculate the ES_Matrix(0_nocom,1_p2p,2_cocom):" matrix_way
read -p "Please enter the way to execute clustering(0_KMediods,1_KMediods++):" cluster_way

# execute matlab script to parse the data  
matlab -nodesktop -nosplash -nojvm -r "file_input='../data/modzs_n272x978.gctx'; file_name='../data/data_for_test.txt'; file_name_cid='../data/data_for_test_cid.txt'; file_name_rid='../data/data_for_test_rid.txt'; PreGSEA; quit;"

outputfile="../data/Cluster_result_test.txt"
n=3
ppn=3
thread_num=5
siglen=50
cluster_num=8

#calculate the similarity(ES) matrix
case "$matrix_way" in
	0)mpirun -n $n -ppn $ppn -hostfile hostfile ES_Matrix_ompi_nocom $thread_num $siglen "../data/data_for_test.txt" "../data/data_for_test.txt" "../data/ES_Matrix_tmp";;
	1)mpirun -n $n -ppn $ppn -hostfile hostfile ES_Matrix_ompi_p2p $thread_num $siglen "../data/data_for_test.txt" "../data/data_for_test.txt" "../data/ES_Matrix_tmp";;
	2)mpirun -n $n -ppn $ppn -hostfile hostfile ES_Matrix_ompi_cocom $thread_num $siglen "../data/data_for_test.txt" "../data/data_for_test.txt" "../data/ES_Matrix_tmp";;
esac

#cluster
case "$cluster_way" in
	0) mpirun -n $n -ppn $ppn -hostfile hostfile Cluster_KMediods_ompi $thread_num $cluster_num "../data/ES_Matrix_tmp" $outputfile;;
	1)mpirun -n $n -ppn $ppn -hostfile hostfile Cluster_KMediods++_ompi $thread_num $cluster_num "../data/ES_Matrix_tmp" $outputfile;;
esac