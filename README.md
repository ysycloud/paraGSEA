# paraGSEA
Gene Set Enrichment Analysis (GSEA) Tools for LINCS data in a parallel manner

## Description

paraGSEA implements MPI and OpenMP-Based parallel GSEA algorithm for multi-core or cluster architecture. 

Many studies have been conducted using gene expression profile similarity to identify functional connections among genes and drugs. While working well for query single gene set in reference of small datasets, its scalability and computation performance is poor in large scale datasets. Here we propose paraGSEA, a parallel computing framework for large scale transcriptomics data analysis. In the process of pairwise similarity metric generation and expression profile clustering, time and space complexity are greatly reduced through elimination of redundant GSEA calculations and optimization of parallel procedures. The framework can be executed on multi-CPU or multi-core computing systems in an efficient manner.

In general, We used the 1ktools(https://github.com/cmap/l1ktools) tools to parse the .gctx file which stored gene profile data defined by Lincs(CMap) based on HDF5 file format.We used Matlab to parse the .gctx file、extract the gene profile sets and write to .txt file. C will read the file 、complete parallel GSEA and write out the result in a suitable file.

## Implementation Details
Several optimizations are implemented in paraGSEA for every stages.

1. In our work, we first implemented GSEA approach in efficient parallel strategy with MPI and OpenMP.In this part, on the one hand, we reduced the computational overhead of standard procedure to calculate the Enrichment Score by pre-sorting, indexing and removing the prefix sum. On the other hand, we will take a global permutation method to wipe off the redundant overhead of estimation of significance level step and multiple hypothesis testing step.

2. Second, we expanded GSEA’s application to quickly compare two gene profile sets to get an Enrichment Score matrix of every gene profile pairs. In this part, in addition to using the previous optimization strategies, our implementation also allows to generate a second level of parallelization by creating several threads per MPI process. The assignment of tasks to threads or processes is performed through a strict load balancing strategy, which leads to a better performance.

3. Third, we clustered the gene profile based on the Enrichment Score matrix which we can get by the second part. In this part, Enrichment Score is served as the metric to measure the similarity between two gene profiles. We implemented a general clustering algorithm like K-Mediods which is an improved version of K-Means. The algorithm can quickly converge and then output the corresponding results.

## Benchmark

With all these optimizations, paraGSEA can attains a 50x speedup compared with original GSEA algorithm in calculating single Enrichment Score. Also, we adopted an global perturbation and random sampling strategy to manage computing expanses in calculate statistical metric of GSEA so that we improved the performance around 100 fold. Moreover, Because of the good data partitioning and communication strategy, the Tools obtained excellent scalability. If the amount of data is large enough, The Tools will keep near linear speedup as the increase of computing nodes.

## Install

### Prerequisites
* [1ktools](https://github.com/cmap/l1ktools) is used to parse the .gctx file which stored gene profile data defined by LINCS and Connectivity Map (CMap) in HDF5 file format.
* Matlab to parse the .gctx file、extract the gene profile sets and export to plain text file, which will be taken as input of C code of parallel GSEA.
* MPICH2
* GCC compiler supports the OpenMP v2.5, v4.0 specification.

### Download and setup

1. You can use the command line `git clone https://github.com/ysycloud/paraGSEA.git` to download the software easily.

2. Configure `Matlab Tools`: you may need to set `paraGSEA/matlab_for_parse` as the `Matlab path` to parse the original file first.

3. Install and configure `C Tools`: In order to install the `C Tools`, you must enter the `paraGSEA` directory first. Then, you can compile and install the Tools in current `bin` directory using command line `make all` so that you can use these tools in this directory. if you want to use the tools in every places of this system, you should run `make install`. However, you may need root authority to execute this operation. And, if you finshed this operation, you can run `make clean` to clean the local Tools in `bin` directory.

### Test

You can test any Tools we provided with some simple datas in `data` directory easily. For example , run `quick_search_serial "data/data_for_test.txt" 10` to test the `Serial Quick Search Tools` in single node. Also, Some typical test sample is provided in [runParaGSEALinux.sh](runParaGSEALinux.sh) shell script. 

The whole install and test process is provided in [install.sh](install.sh) shell script

## Description of Files

### Matlab Tools: matlab_for_parse/

#### Requirements:

1. Matlab R2009a and above

#### Setting the MATLAB path:
Run `cd paraGSEA/matlab_for_parse` or enter the `pathtool` command, click "Add with Subfolders...", and select the directory `paraGSEA/matlab_for_parse`.

#### using by command line:
After you set the MATLAB path, you should enter `matlab` in shell to start matlab environment.
Then, you can set the input and output file path in matlab environment and enter `PreGSEA` to parse original data.
Or, you can use the shell script below to start matlab environment and parse original data.

```shell
% execute matlab script to parse the data  
matlab -nodesktop -nosplash -nojvm -r "file_input='../data/modzs_n272x978.gctx'; file_name='../data/data_for_test.txt'; file_name_cid='../data/data_for_test_cid.txt'; file_name_rid='../data/data_for_test_rid.txt'; PreGSEA; quit;"
```

**Note:** the example of setting input and output file path in matlab script and parse original data is shown below.
```matlab
% setting input file name（ .gctx ----original gene dataset）
file_input = '../data/modzs_n272x978.gctx';

% setting the out file name of profiles
file_name = '../data/data_for_test.txt';  

% setting out file name of profiles' cids(profiles)
file_name_cid = '../data/data_for_test_cid.txt';

% setting out file name of profiles' rids(genes)
file_name_rid='../data/data_for_test_rid.txt';   

% excuting the PreAnalysis in GSEA for C Tools
PreGSEA
```

By the way, there are another more efficent script provided to parse the original data in a parallel manner.
To use this script, you are supposed to make sure that you have a multicores system first, and except the input and output file path you must set in matlab environment like the last script, the number of cores are also should be setted. Then you
can enter `paraPreGSEA` to parse original data.
Or, you can use the shell script below to start matlab environment and parse original data in a parallel manner.

```shell
# execute matlab script to parallel parse the data  
matlab -nodesktop -nosplash -nojvm -r "file_input='../data/modzs_n272x978.gctx'; file_name='../data/data_for_test.txt'; file_name_cid='../data/data_for_test_cid.txt'; file_name_rid='../data/data_for_test_rid.txt'; cores = 2; paraPreGSEA; quit;"

cat ../data/data_for_test.txt_* >> ../data/data_for_test.txt
cat ../data/data_for_test_cid.txt_* >> ../data/data_for_test_cid.txt
rm -f ../data/data_for_test.txt_* ../data/data_for_test_cid.txt_*
```

**Note:** the number of cores must be smaller than the actual core number in your system. And after the parse work, you shoul merge every parts of output file into a whole file like the shell script shown above.

#### Tools:
* [**PreGSEA.m**](matlab_for_parse/PreGSEA.m) : extract the gene profile sets、finish pre-sorting and write to .txt file.
* [**paraPreGSEA.m**] (matlab_for_parse/paraPreGSEA.m): extract the gene profile sets、finish pre-sorting and write to .txt file in a parallel manner.
* [**parse_gctx.m**] (matlab_for_parse/parse_gctx.m): parse .gctx file which is provided by 1ktools.

#### Note:
 * the example of executing shell script to parse the data is provided by [example/runPreGSEAbyMatlab.sh](docs/example/runPreGSEAbyMatlab.md)
 * the example of executing shell script to parse the data in a parallel manner is provided by [example/runparaPreGSEAbyMatlab.sh](docs/example/runparaPreGSEAbyMatlab.md)

### C Tools: src/

#### Requirements:

1. MPICH2
2. GCC compiler supports the OpenMP v2.5, v4.0 specification

#### INSTALL:
* [**install.sh**](install.sh): shell script for Installing all C tools. However, you may need root authority to execute the whole script.

Or, you can use the shell script below easily. Also, you need root authority to run `make install`.
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

* [**quick_search_serial.c**](docs/Tools/quick_search_serial.md) read the .txt file 、complete GSEA and show the topN results in a serial way.
* [**quick_search_omp.c**](docs/Tools/quick_search_omp.md) read the .txt file 、complete parallel GSEA by OpenMP and show the topN results.
* [**quick_search_mpi.c**](docs/Tools/quick_search_mpi.md) read the .txt file 、complete parallel GSEA by MPI and show the topN results.
* [**ES_Matrix_ompi_nocom.c**](docs/Tools/ES_Matrix_ompi_nocom.md) read the .txt file 、complete parallel computing ES_Matrix and write out the result by MPI/OpenMP with no communication in distributing second file of input.
* [**ES_Matrix_ompi_p2p.c**](docs/Tools/ES_Matrix_ompi_p2p.md) read the .txt file 、complete parallel computing ES_Matrix and write out the result by MPI/OpenMP with p2p communication in distributing second file of input.
* [**ES_Matrix_ompi_cocom.c**](docs/Tools/ES_Matrix_ompi_cocom.md) read the .txt file 、complete parallel computing ES_Matrix and write out the result by MPI/OpenMP with collective communication in distributing second file of input.
* [**Cluster_KMediods_ompi.c**](docs/Tools/Cluster_KMediods_ompi.md) read the ES_Matrix file 、complete a general clustering algorithm like K-Mediods by MPI/OpenMP.
* [**Cluster_KMediods++_ompi.c**](docs/Tools/Cluster_KMediods++_ompi.md) read the ES_Matrix file 、complete a general clustering algorithm like K-Mediods by MPI/OpenMP, but let the distance between initial cluster centers far as possible.

#### Demo:
* [**runParalESLinux.sh**](runParaGSEALinux.sh): run Executable files in Linux.

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
 * `example/cluster_demo.sh` executes a whole cluster process to avoid input file and MPI Settings mismatched error which will be mentioned in `Using Problem` part.


### Example
* [**runPreGSEAbyMatlab.sh**](docs/example/runPreGSEAbyMatlab.md): a shell script example to execute pretreatment of parse the original .gctx.
* [**runparaPreGSEAbyMatlab.sh**](docs/example/runparaPreGSEAbyMatlab.md): a shell script example to execute pretreatment of parse the original .gctx in a parallel manner.
* [**quick_search_demo.sh**](docs/example/quick_search_demo.md): a shell script example to execute a whole quick_search process includes parses original data
, select quick search way and quick search.
* [**cluster_demo.sh**](docs/example/cluster_demo.md): a shell script example to execute a whole cluster process includes parses original data
, select ES_Matrix & cluster way and execute ES_Matrix & cluster.

### Datasets in examples

| File | Description | Format|
| ---- | ----------- |--------|
| modzs_n272x978.gctx | original profile file from LINCS Dataset| HDF5 |
| data_for_test.txt | ranked profile file in a index format | first line : profile_number & profile_Length ; next profile_number lines : a ranked profile file included profile_Length elements |
| data_for_test_cid.txt | cid( profile identification ) file | each line : a cid( profile identification ) string corresponding to the last profile_number lines of the `data_for_test.txt` |
| data_for_test_rid.txt | rid( gene identification ) file | each line : a rid( gene identification ) string corresponding to the rid attribute of `modzs_n272x978.gctx` and its index corresponding to the ranked profile file |
| ES_Matrix_test_*.txt | ES Matrix file stored in distributed way ( ‘*’ will be replaced by process id )| first line : row_number	column_number ; next row_number lines : a Enrichment scores vector included column_number elements |
| Cluster_result_test.txt | cluster flag vector file | each line : a cluster flag corresponding to each profile |

## Notes
1. Because of the inefficient IO of Matlab, when the original profile file(.gctx) is too large, the pretreatment operation may take a long time, and it does not support parallel. You may need to be patient. Moreover, Once parsed, it can be reused many times.
2. the program needs a GeneSet as an input in quick_search part after loading the file. You should input a integer string split by space. Each integer in this string represents a gene which you can query in `data_for_test_rid.txt`.
3. When we want to execute the Cluster operator, we must note that input matrix should include the same identity of rows and columns, which means the program that calculates ES Matrix is supposed to use same two file as its input. Only in this way can we get the similarity of each profile pair.
4. When we want to execute the Cluster operator, we must also note that the MPI Settings and hostfile should not be changed compared to the program that calculates ES Matrix. Because the ES_Matrix is stored in distributed way, if you change these settings, each process can not find the right ES matrix blocks. Therefore, if you want to avoid problem 2 and problem 3, you can easily execute the `example/cluster_demo.sh`.
5. If you set the number of clusters too big, clustering algorithm may not converge quickly.

## Connect
Any Question could be sent to the following E-mails:

pittacus@gmail.com, pengshaoliang1979@163.com, cloudysy109@126.com

