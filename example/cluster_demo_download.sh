#!/bin/bash

read -p "Please enter total number of processes:" n
read -p "Please enter the number of processes in per node:" ppn
read -p "Please enter the number of threads in per process:" thread_num
read -p "Please enter the length of gene signature:" siglen
read -p "Please enter the number of clusters:" cluster_num
read -p "Please enter the input profiles file（default by no enter）:" inputfile
read -p "Please enter the output cluster flag vector file:" outputfile
read -p "Please enter the way to calculate the ES_Matrix(0_nocom,1_p2p,2_cocom):" matrix_way
read -p "Please enter the way to excute clustering(0_KMediods,1_KMediods++):" cluster_way

if["$inputfile"=""];then
exit
fi


case "$matrix_way" in
	0) mpirun -n $n -ppn $ppn -hostfile hostfile ES_Matrix_ompi_nocom $thread_num $siglen $inputfile $inputfile "data/ES_Matrix_tmp";;
	1)mpirun -n $n -ppn $ppn -hostfile hostfile ES_Matrix_ompi_p2p $thread_num $siglen $inputfile $inputfile "data/ES_Matrix_tmp";;
	2)mpirun -n $n -ppn $ppn -hostfile hostfile ES_Matrix_ompi_cocom $thread_num $siglen $inputfile $inputfile "data/ES_Matrix_tmp";;
esac

case "$cluster_way" in
	0) mpirun -n $n -ppn $ppn -hostfile hostfile Cluster_KMediods_ompi $thread_num $cluster_num "data/ES_Matrix_tmp" $outputfile;;
	1)mpirun -n $n -ppn $ppn -hostfile hostfile Cluster_KMediods++_ompi $thread_num $cluster_num "data/ES_Matrix_tmp" $outputfile;;
esac