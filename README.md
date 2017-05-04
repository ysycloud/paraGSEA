# paraGSEA
Gene Set Enrichment Analysis (GSEA) Tools for LINCS data in a parallel manner

## Contents: ##
    1. Description.
    2. Benchmark.
    3. Compilation and installation.
    4. Description of Tools.
		4.1 Matlab Tools: matlab_for_parse/
		4.2 C Tools: src/
    5. Datasets in examples | File formats.
	6. Using Problems
    7. License.
    8. Contact.

## I. Description

paraGSEA implements a MPI and OpenMP-Based parallel GSEA algorithm for multi-core or cluster architecture. 

Many studies have been conducted using gene expression profile similarity to identify functional connections among genes and drugs. While working well for query single gene set in reference of small datasets, its scalability and computation performance is poor in large scale datasets. Here we propose paraGSEA, a parallel computing framework for large scale transcriptomics data analysis. In the process of pairwise similarity metric generation and expression profile clustering, time and space complexity are greatly reduced through elimination of redundant GSEA calculations and optimization of parallel procedures. The framework can be executed on multi-CPU or multi-core computing systems in an efficient manner.

In general, We used the 1ktools(https://github.com/cmap/l1ktools) tools to parse the .gctx file which stored gene profile data defined by Lincs(CMap) based on HDF5 file format. We used Matlab to parse the .gctx(.gct) file、extract the gene profile sets and write to .txt file. C will read the file 、complete parallel GSEA and write out the results in some suitable files.

There mainly parts of work and several optimizations are implemented in paraGSEA.

1. First, we implement GSEA approach in efficient parallel strategy with MPI and OpenMP to perform a quick search task, which needs users input a gene set and it will output the top N results after searching the profile data set by carrying out GSEA calculations. In this part, on the one hand, we reduced the computational overhead of standard procedure to calculate the Enrichment Score by pre-sorting, indexing and removing the prefix sum. On the other hand, we will take a global permutation method to wipe off the redundant overhead of estimation of significance level step and multiple hypothesis testing step.

2. Second, we expanded GSEA’s application to quickly compare two gene profile sets to get an Enrichment Score matrix of every gene profile pairs. In this part, in addition to using the previous optimization strategies, our implementation also allows to generate a second level of parallelization by creating several threads per MPI process. The assignment of tasks to threads or processes is performed through a strict load balancing strategy, which leads to a better performance.

3. Third, we clustered the gene profile based on the Enrichment Score matrix which we can get by the second part. In this part, Enrichment Score is served as the metric to measure the similarity between two gene profiles. We implemented a general clustering algorithm like K-Mediods which is an improved version of K-Means. The algorithm can quickly converge and then output the corresponding results. Also, we improved algorithm and provided an implementation of k-mediods++, which is able to ensure that the mutual distances between initial clustering centers as far as possible to achieve better results.

## II. Benchmark

With all these optimizations, paraGSEA can attains a 50x speedup compared with original GSEA algorithm in calculating single Enrichment Score. Also, we adopted an global perturbation and random sampling strategy to manage computing expanses in calculate statistical metric of GSEA so that we improved the performance around 100 fold. Moreover, Because of the good data partitioning and communication strategy, the Tools obtained excellent scalability. If the amount of data is large enough, The Tools will keep near linear speedup as the increase of computing nodes.

## III. Compilation and installation

### III.I. Prerequisites
* [1ktools](https://github.com/cmap/l1ktools) is used to parse the .gctx or .gct file which stored gene profile data defined by LINCS and Connectivity Map (CMap) in HDF5 file format.
* Matlab to parse the .gctx(.gct) file、generate some reference data, extract the gene profile sets and write out to plain text files, which will be taken as input of C code of parallel GSEA.
* MPICH2
* GCC compiler supports the OpenMP v2.5, v4.0 specification.

### III.II. Download and setup

1. You can use the command line `git clone https://github.com/ysycloud/paraGSEA.git` to download the software easily.

2. Configure `Matlab Tools`: you may need to set `paraGSEA/matlab_for_parse` as the `Matlab path` to parse the original file first.

3. Install and configure `C Tools`: In order to install the `C Tools`, you must enter the `paraGSEA` directory first. Then, you can compile and install the Tools in current `bin` directory using command line `make all` so that you can use these tools in this directory. if you want to use the tools in every places of this system, you should run `make install`. However, you may need root authority to execute this operation. Also, if you finshed this operation, you can run `make clean` to clean the local Tools in `bin` directory.

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

### III.III. Test

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


## V. Description of Tools

### V.I. Matlab Tools: matlab_for_parse/

#### V.I.I. Requirements:

1. Matlab R2009a and above

#### V.I.II. Setting the MATLAB path:
Run `cd paraGSEA/matlab_for_parse` or enter the `pathtool` command, click "Add with Subfolders...", and select the directory `paraGSEA/matlab_for_parse`.
Or, if you cannot use Matlab by a visual way, you can just run `addpath('paraGSEA/matlab_for_parse')` after setup the Matlab environment to add the directory path.

#### V.I.III. using by command line:
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

#### V.I.V. Tools:

* [**genReferenceforNewDataSet.m**](matlab_for_parse/genReferenceforNewDataSet.m) : generate some reference data for new data set to facilitate main work of C Tools .
* [**PreGSEA.m**](matlab_for_parse/PreGSEA.m) : extract the gene profile sets、finish pre-sorting and write to .txt file.
* [**paraPreGSEA.m**] (matlab_for_parse/paraPreGSEA.m): extract the gene profile sets、finish pre-sorting and write to .txt file in a parallel manner.
* [**parse_gctx.m**] (matlab_for_parse/lib/parse_gctx.m): parse .gctx file which is provided by 1ktools.
* [**parse_gct.m**] (matlab_for_parse/lib/parse_gct.m): parse .gctx file which is provided by 1ktools.

#### V.I.VI. Demo:
 * the example of executing shell script to generate reference data is provided by [example/runGenReference.sh](docs/example/runGenReference.md)
 * the example of executing shell script to parse the data is provided by [example/runPreGSEAbyMatlab.sh](docs/example/runPreGSEAbyMatlab.md)
 * the example of executing shell script to parse the data in a parallel manner is provided by [example/runparaPreGSEAbyMatlab.sh](docs/example/runparaPreGSEAbyMatlab.md)

### V.II. C Tools: src/

#### V.II.I. Requirements:

1. MPICH2
2. GCC compiler supports the OpenMP v2.5, v4.0 specification

#### V.II.II. Tools(Source file list):

* [**quick_search_serial.c**](docs/Tools/quick_search_serial.md) read the .txt file 、complete GSEA and show the topN results in a serial way.
* [**quick_search_omp.c**](docs/Tools/quick_search_omp.md) read the .txt file 、complete parallel GSEA by OpenMP and show the topN results.
* [**quick_search_mpi.c**](docs/Tools/quick_search_mpi.md) read the .txt file 、complete parallel GSEA by MPI and show the topN results.
* [**ES_Matrix_ompi_nocom.c**](docs/Tools/ES_Matrix_ompi_nocom.md) read the .txt file 、complete parallel computing ES_Matrix and write out the results by MPI/OpenMP with no communication.
* [**ES_Matrix_ompi_p2p.c**](docs/Tools/ES_Matrix_ompi_p2p.md) read the .txt file 、complete parallel computing ES_Matrix and write out the results by MPI/OpenMP with point to point communication.
* [**ES_Matrix_ompi_cocom.c**](docs/Tools/ES_Matrix_ompi_cocom.md) read the .txt file 、complete parallel computing ES_Matrix and write out the results by MPI/OpenMP with collective communication.
* [**Cluster_KMediods_ompi.c**](docs/Tools/Cluster_KMediods_ompi.md) read the ES_Matrix file 、complete a general clustering algorithm like K-Mediods by MPI/OpenMP.
* [**Cluster_KMediods++_ompi.c**](docs/Tools/Cluster_KMediods++_ompi.md) read the ES_Matrix file 、complete a general clustering algorithm like K-Mediods by MPI/OpenMP, but let the distance between initial cluster centers as far as possible.

#### V.II.III. Demo:
* [**runParalESLinux.sh**](runParaGSEALinux.sh): run Executable files in Linux.

#### V.II.V. Running paraGSEA:

paraGSEA runs on Linux and Mac.

List of arguments:

  > quick_search_* [options] [-i INPUT_FILE] [-n TOP_N_RECORDS] [-t THREAD_NUMBER] [-s SAMPLE_SEQUENCE_NUMBER_FILE] [-r REFERENCE_DATA_DIRECTORY]
  
  **-i or --input** *input file*

     a parsed profiles's file from pretreatment stage as input.

  **-t or --thread** *threads*

     Defines the maximum number of parallel threads.
	 must be a positive value

  **-n or --topn** *top N*

	 Define the first and last N GSEA records ordered by ES.
	 must be a positive value
	 
  **-s or --sample** *sample sequece number file*

	 a text file include sample sequence numbers which are extracted from pretreatment stage.

  **-r or --reference** *reference data directory*

	 a directory include some reference data files we generate from pretreatment stage
	 

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


  > mpirun [options] [-n PROCESS_NUM] [-ppn PERNUM] [-hostfile HOSTFILE] Cluster_KMediods*_ompi [options] [-i INPUT_FILE] [-t THREAD_NUMBER] [-c CLUSTER_NUMBERS] [-o OUTPUT_FILE] [-s SAMPLE_SEQUENCE_NUMBER_FILE] [-r REFERENCE_DATA_DIRECTORY]
  
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
	 
  **-s or --sample** *sample sequece number file*

	 a text file include sample sequence numbers which are extracted from pretreatment stage.

  **-r or --reference** *reference data directory*

	 a directory include some reference data files we generate from pretreatment stage

	 
the detail usage of each C Tool is shown below.
```shell
#param list :filename topn
#quick_search_serial -i data/data_for_test.txt -n 10 -s data/data_for_test_cidnum.txt -r data/Reference

#param list :filename thread_num topn
#quick_search_omp -i data/data_for_test.txt -t 4 -n 10 -s data/data_for_test_cidnum.txt -r data/Reference 

#param list :process_num pernum hostfile filename topn
#mpirun -n 2 -ppn 2 -hostfile example/hostfile quick_search_mpi -i data/data_for_test.txt -n 15 -s data/data_for_test_cidnum.txt -r data/Reference

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile ES_Matrix_ompi_nocom -t 4 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix_test

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile ES_Matrix_ompi_p2p -t 4 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix_test

#param list :process_num pernum hostfile thread_num siglen filename1 filename2 outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile ES_Matrix_ompi_cocom -t 4 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix_test

#param list :process_num pernum hostfile thread_num cluster_num filename outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile Cluster_KMediods_ompi -t 4 -c 10 -i data/ES_Matrix_test -o data/Cluster_result_test.txt -s data/data_for_test_cidnum.txt -r data/Reference

#param list :process_num pernum hostfile thread_num cluster_num filename outfilename
#mpirun -n 2 -ppn 2 -hostfile example/hostfile Cluster_KMediods++_ompi -t 4 -c 10 -i data/ES_Matrix_test -o data/Cluster_result_test.txt -s data/data_for_test_cidnum.txt -r data/Reference
```

**Note:** 
 * runParalESLinux.sh annotate a list of execution case of C tools. Removing the annotation, you can using it easily.
 * the details of parameter list of each C tool can be seen in runParalESLinux.sh or the `tutorial`.

#### V.II.VI. Demo
* [**quick_search_demo.sh**](docs/example/quick_search_demo.md): a shell script example to execute a whole quick_search process includes parses original data
, select quick search way and quick search.
* [**cluster_demo.sh**](docs/example/cluster_demo.md): a shell script example to execute a whole cluster process includes parses original data
, select ES_Matrix & cluster way and execute ES_Matrix & cluster.


## VI. Datasets in examples | File Formats

| File | Description | Format|
| ---- | ----------- |--------|
| modzs_n272x978.gctx | original profile file from LINCS Dataset| HDF5 |
| Gene_List.txt | all gene names of every profile in original order recorded in HDF5 source file | one gene name(symbol) per line |
| Samples_Condition.txt | treatment conditions of all profiles in original order recorded in HDF5 source file | one profile's conditions per line |
| Samples_RowByteOffset.txt | Bytes offsets of every line in `Samples_Condition.txt` | every offset value is splitted by `\t` |
| data_for_test.txt | ranked profiles file in a sequence number format | first line : profile_number & profile_Length ; next profile_number lines : a ranked profile includes profile_Length genes in a sequence number format |
| data_for_test_cidnum.txt | profile sequence number file corresponding to these profiles we extract in `data_for_test.txt` | one sequence number per line |
| GeneSet.txt | GeneSet example file | one gene name(symbol) per line |
| ES_Matrix_test_*.txt | ES Matrix file stored in distributed way ( ‘*’ will be replaced by process id )| first line : row_number & column_number ; next row_number lines : a Enrichment scores vector included column_number elements |
| Cluster_result_test.txt | cluster result file | each line consists of a class label or a profile information |


### VI.I generated reference data format:
  
  When we get a new source file in correct HDF5 format to analysis, such as the example `modzs_n272x978.gctx`, 
  we need generate some reference data first. 
  Three files will generate as reference data(`Gene_List.txt`, `Samples_Condition.txt`, `Samples_RowByteOffset.txt`).
  
####  VI.I.I. `Gene_List.txt`: ####
  
  This file includes all gene names of every profile in original order recorded in HDF5 source file.
  When users input a GeneSet, we can get the sequence number of every gene with this file.
  The main format is one gene name(symbol) per line.
  
  Example:
	
	PSME1
	ATF1
	RHEB
	FOXO3
	RHOA
	IL1B
	ASAH1
	......
	
####  VI.I.II. `Samples_Condition.txt`: ####
  
  This file includes treatment conditions of all profiles in original order recorded in HDF5 source file.
  When we get a profile sequence number, we can get the detail information of this profile treatment conditions with this file.
  The main format is one profile's conditions per line.

  Example:
  
	cid:CPC006_A549_6H:BRD-U88459701-000-01-8:10;    cell_line:      A549;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	cid:CPC020_A375_6H:BRD-A82307304-001-01-8:10;    cell_line:      A375;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	cid:CPC020_HT29_6H:BRD-A82307304-001-01-8:10;    cell_line:      HT29;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	......

####  VI.I.III. `Samples_RowByteOffset.txt`: ####
  
  This file includes Bytes offsets of every line in `Samples_Condition.txt`.
  With this file, we can locate specific line treatment conditions directly without loading all `Samples_Condition.txt` into memory, 
  that will be very effective in time and space consumption.
  The main format is every offset value is splitted by `\t`.

  Example:
	
	  0	       189	       378	       567	       756	       945	      1135	      1325	......
  

### VI.II Standard Input format of Quick Search:

####  VI.II.I. `data_for_test.txt`: ####
  
  This file is ranked profiles file which is parsed and extracted from source HDF5 file.
  There are two parts of this file. The first is the first line, which only consists of two figures, the profile number and profile length.
  The second part is the content of every profile. Each line is a profile includes profile_Length genes in a sequence number format.
  In order to get this part, we first numbered every gene starting with one and then ordered them according to their differential expression.
 
  Example:
  
	272	       978
	562	  832	  688	  690	  136	  895	  682	...
	650	  605	  711	  436	  630	  429	  787	...  
	175	  136	  857	  145	  832	  707	  850	...
	 80	  207	  102	  127	  861	  512	  860	...
	......

####  VI.II.II. `data_for_test_cidnum.txt`: ####
  
  This file is profile sequence number file corresponding to these profiles we extract in `data_for_test.txt`.
  When we get this profile sequence number, we can get the detail information of this profile treatment conditions with reference data.
  The main format is one sequence number per line.
 
  Example:
  
        14
        15
        18
        19
        20
        21
		......
		
####  VI.II.III. `GeneSet.txt`: ####
	
  This is a GeneSet example file with the same format of `Gene_List.txt`, where presents as one gene name(symbol) per line.
	
  Example:
  
	CDKN2A
	CDKN1B
	GAPDH
	CISD1
	SPDEF
	IGF1R
	......
	
### VI.III Standard Input format of Compare Profiles:
  
  The mainly input file of this part is two ranked profiles files, such as `data_for_test.txt`, 
  which has been described before. Therefore, there will be no more description.
	

### VI.V Standard Input format of Clusting Profiles( Output format of Compare Profiles ):

####  VI.V.I. `ES_Matrix_test_*.txt`: ####

  There are some ES Matrix files stored in distributed way ( ‘*’ will be replaced by process id ). 
  There are also two parts of each ES Matrix file. The first is the first line, which only consists of two figures, the row number and column number.
  The second part is the ES vectors of every profile to all other profiles, where each line is a enrichment scores vector included column_number elements.
  
  Example:
  
	136	       272
	1.000	0.271	0.147	0.067	0.247	-0.065 ...
	0.271	1.000	0.259	0.109	0.265	0.256  ...
	0.147	0.259	1.000	0.071	-0.012	-0.061 ...
	0.067	0.109	0.071	1.000	0.185	0.433  ...
	0.247	0.265	-0.012	0.185	1.000	0.226  ...
	......

### VI.VI Standard Output format of Clusting Profiles:

####  VI.VI.I. `Cluster_result_test.txt`: ####

  This is the cluster result file, which includes the class labels and corresponding profile treatment conditions information.
  
  Example:

	cluster 1 :
	cid:CPC006_A549_6H:BRD-U88459701-000-01-8:10;    cell_line:      A549;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	cid:CPC020_A375_6H:BRD-A82307304-001-01-8:10;    cell_line:      A375;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	......
	cluster 2 :
	cid:CPC020_HT29_6H:BRD-A82307304-001-01-8:10;    cell_line:      HT29;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	cid:CPC006_A375_6H:BRD-U88459701-000-01-8:10;    cell_line:      A375;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	......
	cluster 3 :
	cid:CPC006_A549_24H:BRD-U88459701-000-01-8:10;    cell_line:      A549;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:      24 h;    concentration:     10 uM
	cid:CPC006_AGS_6H:BRD-U88459701-000-01-8:10;    cell_line:       AGS;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	......
	cluster 4 :
	cid:CPC006_HA1E_6H:BRD-U88459701-000-01-8:10;    cell_line:      HA1E;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	cid:CPC006_HCC15_6H:BRD-U88459701-000-01-8:10;    cell_line:     HCC15;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	......
	......
	

## VII. Using Problem
1. When we get a new profile file keeps in correct format with a ‘gctx’ or ‘gct’ suffix to analysis, we must generate some reference data first to facilitate the main follow-up work of C Tools.
2. Because of the inefficient IO of Matlab, when the original profile file is too large, the pretreatment operation may take a long time. You may need to be patient. However, Once parsed, it can be reused many times. Also, there is parallel way are provided to accelerate the pretreatment operation if the multi-cores environment is supported in your machine.
3. the program needs a GeneSet as an input in quick_search part after loading the file. We provided two ways to support the GeneSet input. One is inputting the GeneSet directly which is splitted by space and the other is inputting a file path where including a GeneSet. Most of time, the second way may be more convenient. 
4. When we want to execute the Cluster operator, we must note that input matrix should include the same identity of rows and columns, which means the program that calculates ES Matrix is supposed to use same two file as its input. Only in this way can we get the similarity of each profile pair.
5. When we want to execute the Cluster operator, we must also note that the MPI Settings and hostfile should not be changed compared to the program that calculates ES Matrix. Because the ES_Matrix is stored in distributed way, if you change these settings, each process can not find the right ES matrix blocks. Therefore, if you want to avoid problem 4 and problem 5, you can easily execute the `example/cluster_demo.sh`.
6. If you set the number of clusters too big, clustering algorithm may not converge quickly.


## VIII. License
paraGSEA is licensed under the GNU General Public License, version 3
(GPLv3), for more information read the LICENSE file or refer to:

  http://www.gnu.org/licenses/

## VIIII. Contact
Any Question could be sent to the following E-mails:

pittacus@gmail.com, pengshaoliang1979@163.com, cloudysy109@126.com
