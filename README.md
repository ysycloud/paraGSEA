# paraGSEA
Gene Set Enrichment Analysis (GSEA) Tools for LINCS data in a parallel manner

## I. Description

paraGSEA implements a MPI and OpenMP-Based parallel GSEA algorithm for multi-core or cluster architecture. 

Many studies have been conducted using gene expression profile similarity to identify functional connections among genes and drugs. While working well for query single gene set in reference of small datasets, its scalability and computation performance is poor in large scale datasets. Here we propose paraGSEA, a parallel computing framework for large scale transcriptomics data analysis. In the process of pairwise similarity metric generation and expression profile clustering, time and space complexity are greatly reduced through elimination of redundant GSEA calculations and optimization of parallel procedures. The framework can be executed on multi-CPU or multi-core computing systems in an efficient manner.

In general, We used the 1ktools(https://github.com/cmap/l1ktools) tools to parse the .gctx file which stored gene profile data defined by Lincs(CMap) based on HDF5 file format.We used Matlab to parse the .gctx file、extract the gene profile sets and write to .txt file. C will read the file 、complete parallel GSEA and write out the result in a suitable file.

## II. Implementation Details
There mainly parts of work and several optimizations are implemented in paraGSEA.

1. First, we implement GSEA approach in efficient parallel strategy with MPI and OpenMP to perform a quick search task, which needs users input a gene set and it will output the top N results after searching the profile data set by carrying out GSEA calculations. In this part, on the one hand, we reduced the computational overhead of standard procedure to calculate the Enrichment Score by pre-sorting, indexing and removing the prefix sum. On the other hand, we will take a global permutation method to wipe off the redundant overhead of estimation of significance level step and multiple hypothesis testing step.

2. Second, we expanded GSEA’s application to quickly compare two gene profile sets to get an Enrichment Score matrix of every gene profile pairs. In this part, in addition to using the previous optimization strategies, our implementation also allows to generate a second level of parallelization by creating several threads per MPI process. The assignment of tasks to threads or processes is performed through a strict load balancing strategy, which leads to a better performance.

3. Third, we clustered the gene profile based on the Enrichment Score matrix which we can get by the second part. In this part, Enrichment Score is served as the metric to measure the similarity between two gene profiles. We implemented a general clustering algorithm like K-Mediods which is an improved version of K-Means. The algorithm can quickly converge and then output the corresponding results.

## III. Benchmark

With all these optimizations, paraGSEA can attains a 50x speedup compared with original GSEA algorithm in calculating single Enrichment Score. Also, we adopted an global perturbation and random sampling strategy to manage computing expanses in calculate statistical metric of GSEA so that we improved the performance around 100 fold. Moreover, Because of the good data partitioning and communication strategy, the Tools obtained excellent scalability. If the amount of data is large enough, The Tools will keep near linear speedup as the increase of computing nodes.

## V. Compilation and installation

### V.I. Prerequisites
* [1ktools](https://github.com/cmap/l1ktools) is used to parse the .gctx or .gct file which stored gene profile data defined by LINCS and Connectivity Map (CMap) in HDF5 file format.
* Matlab to parse the .gctx(.gct) file、generate some reference data, extract the gene profile sets and write out to plain text files, which will be taken as input of C code of parallel GSEA.
* MPICH2
* GCC compiler supports the OpenMP v2.5, v4.0 specification.

### V.II. Download and setup

1. You can use the command line `git clone https://github.com/ysycloud/paraGSEA.git` to download the software easily.

2. Configure `Matlab Tools`: you may need to set `paraGSEA/matlab_for_parse` as the `Matlab path` to parse the original file first.

3. Install and configure `C Tools`: In order to install the `C Tools`, you must enter the `paraGSEA` directory first. Then, you can compile and install the Tools in current `bin` directory using command line `make all` so that you can use these tools in this directory. if you want to use the tools in every places of this system, you should run `make install`. However, you may need root authority to execute this operation. Also, if you finshed this operation, you can run `make clean` to clean the local Tools in `bin` directory.

### V.III. Test

You can test any Tools we provided with some simple datas in `data` directory easily. For example , run `quick_search_serial` to test the `Serial Quick Search Tools` in single node, which wiil obtain the output shown below. Also, Some typical test sample is provided in [runParaGSEALinux.sh](runParaGSEALinux.sh) shell script. 
```shell
Usage:   quick_search_serial [options]
	general options:
		-n --topn: The first and last N GSEA records ordered by ES. [ default 10]
	input/output options:
		-i --input: input file/a parsed profiles file from pretreatment stage.
		-s --sample: input file/a parsed sample sequence number file from pretreatment stage.
		-r --reference: input a directory includes referenced files about genesymbols and cids.
```
The whole install and test process is provided in [install.sh](install.sh) shell script

## VI. Description of Tools

### VI.I. Matlab Tools: matlab_for_parse/

#### VI.I.I. Requirements:

1. Matlab R2009a and above

#### VI.I.II. Setting the MATLAB path:
Run `cd paraGSEA/matlab_for_parse` or enter the `pathtool` command, click "Add with Subfolders...", and select the directory `paraGSEA/matlab_for_parse`.

#### VI.I.III. using by command line:
After you set the MATLAB path, you should enter `matlab` in shell to start matlab environment.

In order to provide user-friendly parsed method to allow user set their own conditions of profile they need, we must generate some reference data to facilitate our main work. There is a Matlab script in ‘paraGSEA/matlab_for_parse’ directory named `genReferenceforNewDataSet.m` to help us finish this work. The Only thing we need to do is just setting some field names and file path. There is a example below. More detail have been told in `tutorial`. 

```shell
datasource='../data/modzs_n272x978.gctx';
gene_symbol_rhd = 'pr_gene_symbol'; 
sample_conditions_chd = {'cell_id', 'pert_iname', 'pert_type', 'pert_itime', 'pert_idose'};
genReferenceforNewDataSet
```

Then, you can set the input,output file path and some sample conditions in matlab environment and enter `PreGSEA` to parse original data.
Or, you can use the shell script below to start matlab environment and parse original data directly.

```shell
% execute matlab script to parse the data  
matlab -nodesktop -nosplash -nojvm -r " file_input='../data/modzs_n272x978.gctx';  file_name='../data/data_for_test.txt';  file_name_cidnum='../data/data_for_test_cidnum.txt';  sample_conditions_chd = {'cell_id', 'pert_iname', 'pert_type', 'pert_itime', 'pert_idose'};  cell_id_set={'A549','MCF7','A375','A673','AGS'};  pert_set={'atorvastatin','vemurafenib','venlafaxine'}; pert_type_set = {'trt_cp'}; duration = 6 ;  concentration= 10;  PreGSEA;  quit;"

cat ../data/tmp >> ../data/data_for_test.txt
rm -f ../data/tmp
```

By the way, there are another more efficent script provided to parse the original data in a parallel manner.
To use this script, you are supposed to make sure that you have a multicores system first, and except the input and output file path you must set in matlab environment like the last script, the number of cores are also should be setted. Then you can enter `paraPreGSEA` to parse original data.
Or, you can use the shell script below to start matlab environment and parse original data in a parallel manner.

```shell
# execute matlab script to parallel parse the data  
matlab -nodesktop -nosplash -nojvm -r " file_input='../data/modzs_n272x978.gctx';  file_name='../data/data_for_test.txt';  file_name_cid='../data/data_for_test_cidnum.txt';  cores = 2;  sample_conditions_chd = {'cell_id', 'pert_iname', 'pert_type', 'pert_itime', 'pert_idose'};  cell_id_set={'A549','MCF7','A375','A673','AGS'}; pert_set={'atorvastatin','vemurafenib','venlafaxine'}; pert_type_set = {'trt_cp'}; duration = 6 ;  concentration= 10;  paraPreGSEA;  quit;"

cat ../data/data_for_test.txt_* >> ../data/data_for_test.txt
cat ../data/data_for_test_cidnum.txt_* >> ../data/data_for_test_cidnum.txt
rm -f ../data/data_for_test.txt_* ../data/data_for_test_cidnum.txt_*
```

**Note:** the number of cores must be smaller than the actual core number in your system. And after the parse work, you shoul merge every parts of output file into a whole file like the shell script shown above.

#### VI.I.V. Tools:
* [**PreGSEA.m**](matlab_for_parse/PreGSEA.m) : extract the gene profile sets、finish pre-sorting and write to .txt file.
* [**paraPreGSEA.m**] (matlab_for_parse/paraPreGSEA.m): extract the gene profile sets、finish pre-sorting and write to .txt file in a parallel manner.
* [**parse_gctx.m**] (matlab_for_parse/lib/parse_gctx.m): parse .gctx file which is provided by 1ktools.
* [**parse_gct.m**] (matlab_for_parse/lib/parse_gct.m): parse .gctx file which is provided by 1ktools.

#### VI.I.VI. Note:
 * the example of executing shell script to parse the data is provided by [example/runPreGSEAbyMatlab.sh](docs/example/runPreGSEAbyMatlab.md)
 * the example of executing shell script to parse the data in a parallel manner is provided by [example/runparaPreGSEAbyMatlab.sh](docs/example/runparaPreGSEAbyMatlab.md)

### VI.II. C Tools: src/

#### VI.II.I. Requirements:

1. MPICH2
2. GCC compiler supports the OpenMP v2.5, v4.0 specification

#### VI.II.II. Compilation and installation:
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
```

#### VI.II.III. Tools:

* [**quick_search_serial.c**](docs/Tools/quick_search_serial.md) read the .txt file 、complete GSEA and show the topN results in a serial way.
* [**quick_search_omp.c**](docs/Tools/quick_search_omp.md) read the .txt file 、complete parallel GSEA by OpenMP and show the topN results.
* [**quick_search_mpi.c**](docs/Tools/quick_search_mpi.md) read the .txt file 、complete parallel GSEA by MPI and show the topN results.
* [**ES_Matrix_ompi_nocom.c**](docs/Tools/ES_Matrix_ompi_nocom.md) read the .txt file 、complete parallel computing ES_Matrix and write out the result by MPI/OpenMP with no communication in distributing second file of input.
* [**ES_Matrix_ompi_p2p.c**](docs/Tools/ES_Matrix_ompi_p2p.md) read the .txt file 、complete parallel computing ES_Matrix and write out the result by MPI/OpenMP with p2p communication in distributing second file of input.
* [**ES_Matrix_ompi_cocom.c**](docs/Tools/ES_Matrix_ompi_cocom.md) read the .txt file 、complete parallel computing ES_Matrix and write out the result by MPI/OpenMP with collective communication in distributing second file of input.
* [**Cluster_KMediods_ompi.c**](docs/Tools/Cluster_KMediods_ompi.md) read the ES_Matrix file 、complete a general clustering algorithm like K-Mediods by MPI/OpenMP.
* [**Cluster_KMediods++_ompi.c**](docs/Tools/Cluster_KMediods++_ompi.md) read the ES_Matrix file 、complete a general clustering algorithm like K-Mediods by MPI/OpenMP, but let the distance between initial cluster centers far as possible.

#### VI.II.V. Demo:
* [**runParalESLinux.sh**](runParaGSEALinux.sh): run Executable files in Linux.

#### VI.II.VI. Running paraGSEA:

paraGSEA runs on Linux and Mac.

List of arguments:

  > quick_search_* [options] [-i INPUT_FILE] [-n TOP_N_RECORDS] [-t THREAD_NUMBER]
  
  **-i or --input** *input file*

     a parsed profiles's file from pretreatment stage as input.

  **-t or --thread** *threads*

     Defines the maximum number of parallel threads.
	 must be a positive value

  **-n or --topn** *top N*

	 Define the first and last N GSEA records ordered by ES.
	 must be a positive value
	 

  > mpirun [options] [-n PROCESS_NUM] [-ppn PERNUM] [-hostfile HOSTFILE] ES_Matrix_ompi_* [options] [-1 INPUT_FILE1] [-2 INPUT_FILE2] [-l SIGLEN] [-t THREAD_NUMBER] [-o OUTPUT_FILE]
  
  **-n** *mpi parameter*

      Total number of processes.
	   
  **-ppn** *mpi parameter*

     the number of processes in each node.

  **-hostfile** *mpi parameter*

     list the IP or Hostname of nodes

	 you'd better keep the formula correct : process_num = pernum * number of IP(Hostname) list in hostfile
  
  **-1 or --input1** *input file1*

     a parsed profiles's file from pretreatment stage as input1.
	   
  **-2 or --input2** *input file2*

     a parsed profiles's file from pretreatment stage as input2.

  **-t or --thread** *threads*

     Defines the maximum number of parallel threads.
	 must be a positive value

  **-l or --siglen** *signature length*

	 Define the length of Gene Expression Signature.
	 must be a positive value	 
	 
  **-o or --output** *output file*
  
	 Define the output file ,distributed in every nodes ,with ES Matrix


  > mpirun [options] [-n PROCESS_NUM] [-ppn PERNUM] [-hostfile HOSTFILE] Cluster_KMediods*_ompi [options] [-i INPUT_FILE] [-t THREAD_NUMBER] [-c CLUSTER_NUMBERS] [-o OUTPUT_FILE]
  
  **-n** *mpi parameter*

      Total number of processes.
	   
  **-ppn** *mpi parameter*

     the number of processes in each node.

  **-hostfile** *mpi parameter*

     list the IP or Hostname of nodes

	 you'd better keep the formula correct : process_num = pernum * number of IP(Hostname) list in hostfile
  
  **-i or --input** *input file*

     distributed ES_Matrix file we get from stage 2(Compare Profiles).

  **-t or --thread** *threads*

     Defines the maximum number of parallel threads.
	 must be a positive value

  **-c or --cluster** *cluster number*

	 Define the number of clusters we want to get.
	 must be a positive value	 
	 
  **-o or --output** *output file*
  
	 Define the output cluster result file of every profiles in root node


the detail usage of each C Tools is shown below.
```shell
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
#mpirun -n 2 -ppn 2 -hostfile example/hostfile Cluster_KMediods_ompi -t 4 -c 10 -i data/ES_Matrix_test -o data/Cluster_result_test.txt

#param list :process_num pernum hostfile thread_num cluster_num filename outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile Cluster_KMediods++_ompi -t 4 -c 10 -i data/ES_Matrix_test -o data/Cluster_result_test.txt
```

#### VI.II.VII. Note:
 * runParalESLinux.sh annotate a list of execution case of C tools. Removing the annotation, you can using it easily.
 * the details of parameter list of each C tools can be seen in runParalESLinux.sh or the TUTORIAL.
 * `example/cluster_demo.sh` executes a whole cluster process to avoid input file and MPI Settings mismatched error which will be mentioned in `Using Problem` part.


### VI.III. Example
* [**runPreGSEAbyMatlab.sh**](docs/example/runPreGSEAbyMatlab.md): a shell script example to execute pretreatment of parse the original .gctx.
* [**runparaPreGSEAbyMatlab.sh**](docs/example/runparaPreGSEAbyMatlab.md): a shell script example to execute pretreatment of parse the original .gctx in a parallel manner.
* [**quick_search_demo.sh**](docs/example/quick_search_demo.md): a shell script example to execute a whole quick_search process includes parses original data
, select quick search way and quick search.
* [**cluster_demo.sh**](docs/example/cluster_demo.md): a shell script example to execute a whole cluster process includes parses original data
, select ES_Matrix & cluster way and execute ES_Matrix & cluster.

### VI.V. Datasets in examples | File Formats

| File | Description | Format|
| ---- | ----------- |--------|
| modzs_n272x978.gctx | original profile file from LINCS Dataset| HDF5 |
| data_for_test.txt | ranked profile file in a index format | first line : profile_number & profile_Length ; next profile_number lines : a ranked profile file included profile_Length elements |
| data_for_test_cid.txt | cid( profile identification ) file | each line : a cid( profile identification ) string corresponding to the last profile_number lines of the `data_for_test.txt` |
| data_for_test_rid.txt | rid( gene identification ) file | each line : a rid( gene identification ) string corresponding to the rid attribute of `modzs_n272x978.gctx` and its index corresponding to the ranked profile file |
| ES_Matrix_test_*.txt | ES Matrix file stored in distributed way ( ‘*’ will be replaced by process id )| first line : row_number	column_number ; next row_number lines : a Enrichment scores vector included column_number elements |
| Cluster_result_test.txt | cluster result file | each line : a cluster flag corresponding to each profile |

#### VI.V.I Standard parsed profile format:
  
  Example:
  
	272	       978
	562	  832	  688	  690	  136	  895	  682	...
	650	  605	  711	  436	  630	  429	  787	...  
	175	  136	  857	  145	  832	  707	  850	...
	 80	  207	  102	  127	  861	  512	  860	...


## VII. Notes
1. Because of the inefficient IO of Matlab, when the original profile file(.gctx) is too large, the pretreatment operation may take a long time, and it does not support parallel. You may need to be patient. Moreover, Once parsed, it can be reused many times.
2. the program needs a GeneSet as an input in quick_search part after loading the file. You should input a integer string split by space. Each integer in this string represents a gene which you can query in `data_for_test_rid.txt`.
3. When we want to execute the Cluster operator, we must note that input matrix should include the same identity of rows and columns, which means the program that calculates ES Matrix is supposed to use same two file as its input. Only in this way can we get the similarity of each profile pair.
4. When we want to execute the Cluster operator, we must also note that the MPI Settings and hostfile should not be changed compared to the program that calculates ES Matrix. Because the ES_Matrix is stored in distributed way, if you change these settings, each process can not find the right ES matrix blocks. Therefore, if you want to avoid problem 2 and problem 3, you can easily execute the `example/cluster_demo.sh`.
5. If you set the number of clusters too big, clustering algorithm may not converge quickly.


## VIII. License
paraGSEA is licensed under the GNU General Public License, version 3
(GPLv3), for more information read the LICENSE file or refer to:

  http://www.gnu.org/licenses/

## VIIII. Contact
Any Question could be sent to the following E-mails:

pittacus@gmail.com, pengshaoliang1979@163.com, cloudysy109@126.com

