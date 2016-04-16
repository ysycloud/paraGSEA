# paraGSEA
Gene set Enrichment analysis(GSEA) Tools for Lincs data in a parallel manner

## Analysis Tools

A brief description of the tools is given below. The tools implement the function of parallel compute GSEA in multi-core environment.I used the 1ktools(https://github.com/cmap/l1ktools) tools to parse the .gctx file which stored gene profile data defined by Lincs(CMap) based on HDF5 file format.I used Matlab to parse the .gctx file、extract the gene profile sets and write to .txt file. C will read the file 、complete parallel GSEA and write out the result in a suitable file. 

### General Function:
In our work, we first implemented GSEA approach in efficient parallel strategy with MPI and OpenMP.In this part, on the one hand, we reduced the computational overhead of standard procedure to calculate the Enrichment Score by pre-sorting, indexing and removing the prefix sum. On the other hand, we will take a global permutation method to wipe off the redundant overhead of estimation of significance level step and multiple hypothesis testing step. 

Second, we expanded GSEA’s application to quickly compare two gene profile sets to get an Enrichment Score matrix of every gene profile pairs. In this part, in addition to using the previous optimization strategies, our implementation also allows to generate a second level of parallelization by creating several threads per MPI process.

Third, we clustered the gene profile based on the Enrichment Score matrix which we can get by the second part. In this part, Enrichment Score is served as the metric to measure the similarity between two gene profiles. We implemented a general clustering algorithm like K-Mediods which is an improved version of K-Means. The algorithm can quickly converge and then output the corresponding results.

I will update the tools as they become available.


### Matlab Tools: matlab_for_parse/

#### Requirements:

1. Matlab R2009a and above

#### Setting the MATLAB path:
Run `cd paraGSEA/matlab_for_parse` or enter the `pathtool` command, click "Add with Subfolders...", and select the directory `paraGSEA/matlab_for_parse`.

#### using by command line:
After you set the MATLAB path, you should enter `matlab` in shell to start matlab environment.
Then, you can set the input and output file path in `runPreGSEA.m` and enter `runPreGSEA` to parse the data.
Or, you can use the shell script below to start matlab environment and parse original data after set the 
input and output file path in `runPreGSEA.m`.

```shell
% execute matlab script to parse the data  
matlab -nodesktop -nosplash - nojvm -r "runPreGSEA; quit;"
```
**Note:** the example of setting input and output file path in `runPreGSEA.m` is shown below.
```matlab
% setting input file name（ .gctx ----original gene dataset）
file_input = '../data/modzs_n272x978.gctx'; 

% setting the out file name of profiles
file_name = '../data/data_for_test.txt';  

% setting out file name of profiles' cids 
file_name_cid = '../data/data_for_test_cid.txt';   

% excuting the PreAnalysis in GSEA for C Tools
PreGSEA

```

#### Tools:
* [**runPreGSEA.m**] : Setting Parameters and execute the PreGSEA.
* [**PreGSEA.m**] : extract the gene profile sets、finish pre-sorting and write to .txt file.

### C Tools: src/

#### Requirements:

1. MPI
2. gcc compiler supports OpenMP

#### INSTALL:
* [**install.sh**]: shell script for Installing all C tools.

Or, you can use the shell script below easily.
```shell
#git clone
git clone https://github.com/ysycloud/paraGSEA.git
cd paraGSEA
#make
make all
#install
make install
#clean
make clean
```

#### Tools:

* [**quick_search_serial.c**] read the .txt file 、complete GSEA and show the topN results in a serial way.
* [**quick_search_omp.c**] read the .txt file 、complete parallel GSEA by OpenMP and show the topN results.
* [**quick_search_mpi.c**] read the .txt file 、complete parallel GSEA by MPI and show the topN results.
* [**ES_Matrix_ompi_nocom.c**] read the .txt file 、complete parallel computing ES_Matrix and write out the result by MPI/OpenMP with no communication in distributing second file of input.
* [**ES_Matrix_ompi_p2p.c**] read the .txt file 、complete parallel computing ES_Matrix and write out the result by MPI/OpenMP with p2p communication in distributing second file of input.
* [**ES_Matrix_ompi_cocom.c**] read the .txt file 、complete parallel computing ES_Matrix and write out the result by MPI/OpenMP with collective communication in distributing second file of input.
* [**Cluster_KMediods_ompi.c**] read the ES_Matrix file 、complete a general clustering algorithm like K-Mediods by MPI/OpenMP.
* [**Cluster_KMediods++_ompi.c**] read the ES_Matrix file 、complete a general clustering algorithm like K-Mediods by MPI/OpenMP, but let the distance between initial cluster centers far as possible.

#### Demo:
* [**runParalESLinux.sh**]: run Executable files in Linux.
* [**cluster_demo.sh**]: a shell script to excute a whole cluster process include ES_Matrix and cluster.

the detail usage of each C Tools is shown below.
```shell
#param list :filename topn
quick_search_serial "data/data_for_test.txt" 10

#param list :filename thread_num topn paraway(0/1)
quick_search_omp "data/data_for_test.txt" 4 10 1

#param list :process_num pernum hostfile filename topn
mpirun -n 2 -ppn 2 -hostfile hostfile quick_search_mpi "data/data_for_test.txt" 15

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
mpirun -n 2 -ppn 2 -hostfile hostfile ES_Matrix_ompi_nocom 4 50 "data/data_for_test.txt" "data/data_for_test.txt" "data/ES_Matrix_test"

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
mpirun -n 2 -ppn 2 -hostfile hostfile ES_Matrix_ompi_p2p 4 50 "data/data_for_test.txt" "data/data_for_test.txt" "data/ES_Matrix_test"

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
mpirun -n 2 -ppn 2 -hostfile hostfile ES_Matrix_ompi_cocom 4 50 "data/data_for_test.txt" "data/data_for_test.txt" "data/ES_Matrix_test"

#param list :process_num pernum hostfile thread_num cluster_num filename outfilename
mpirun -n 2 -ppn 2 -hostfile hostfile Cluster_KMediods_ompi 4 12 "data/ES_Matrix_test" "data/Cluster_result_test.txt"

#param list :process_num pernum hostfile thread_num cluster_num filename outfilename
mpirun -n 2 -ppn 2 -hostfile hostfile Cluster_KMediods++_ompi 4 12 "data/ES_Matrix_test" "data/Cluster_result_test.txt"
```

#### Note:
 * runParalESLinux.sh annotate a list of execution case of C tools. Removing the annotation, you can using it easily.
 * the details of parameter list of each C tools can be seen in runParalESLinux.sh or the TUTORIAL.
 * cluster_demo.sh executes a whole cluster process to avoid inputfile and MPI Settings mismatched error, but there are some parameters must be inputed by hand.
 
 the detail of this script is shown below.
```shell
read -p "Please enter total number of processes:" n
read -p "Please enter the number of processes in per node:" ppn
read -p "Please enter the number of threads in per process:" thread_num
read -p "Please enter the length of gene signature:" siglen
read -p "Please enter the number of clusters:" cluster_num
read -p "Please enter the input profiles file:" inputfile
read -p "Please enter the output cluster flag vector file:" outputfile
read -p "Please enter the way to calculate the ES_Matrix(0_nocom,1_p2p,2_cocom):" matrix_way
read -p "Please enter the way to excute clustering(0_KMediods,1_KMediods++):" cluster_way

case "$matrix_way" in
	0) mpirun -n $n -ppn $ppn -hostfile hostfile ES_Matrix_ompi_nocom $thread_num $siglen $inputfile $inputfile "data/ES_Matrix_tmp";;
	1)mpirun -n $n -ppn $ppn -hostfile hostfile ES_Matrix_ompi_p2p $thread_num $siglen $inputfile $inputfile "data/ES_Matrix_tmp";;
	2)mpirun -n $n -ppn $ppn -hostfile hostfile ES_Matrix_ompi_cocom $thread_num $siglen $inputfile $inputfile "data/ES_Matrix_tmp";;
esac

case "$cluster_way" in
	0) mpirun -n $n -ppn $ppn -hostfile hostfile Cluster_KMediods_ompi $thread_num $cluster_num "data/ES_Matrix_tmp" $outputfile;;
	1)mpirun -n $n -ppn $ppn -hostfile hostfile Cluster_KMediods++_ompi $thread_num $cluster_num "data/ES_Matrix_tmp" $outputfile;;
esac
```
 
### Description of Files appeared in examples

| File | Description | Format|
| ---- | ----------- |--------|
| modzs_n272x978.gctx | original profile file from LINCS Dataset| HDF5 |
| data_for_test.txt | ranked profile file | first line : profile_number	profile_Length ; next profile_number lines : a ranked profile file included profile_Length elements |
| data_for_test_cid.txt | cid file | each line : a cid string corresponding to the last profile_number lines first outputfile  |
| ES_Matrix_test_*.txt | ES Matrix file stored in distributed way ( ‘*’ will be replaced by process id )| first line : row_number	column_number ; next row_number lines : a Enrichment scores vector included column_number elements |
| Cluster_result_test.txt | cluster flag vector file | each line : a cluster flag corresponding to each profile  |

## Using Problem
1. Because of the inefficient IO of Matlab, when the original profile file(.gctx) is too large, the pretreatment operation may take a long time, and it does not support parallel. You may need to be patient. 
2. When we want to excute the Cluster operator, we must note that input matrix should include the same identity of rows and columns, which means the program that calculates ES Matrix is supposed to use same two file as its input. Only in this way can we get the similarity of each profile pair.
3. When we want to excute the Cluster operator, we must also note that the MPI Settings and hostfile should not be changed compared to the program that calculates ES Matrix. Because the ES_Matrix is stored in distributed way, if you change these settings, each process can not find the right ES matrix blocks. Therefore, if you want to avoid problem 2 and problem 3, you can easily execute the cluster_demo.sh.
4. If you set the number of clusters too big, clustering algorithm may not converge quickly.