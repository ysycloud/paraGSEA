# paraGSEA
Gene set Enrichment analysis(GSEA) Tools for Lincs data in a parallel manner

## Analysis Tools

A brief description of the tools is given below. The tools implement the function of parallel compute GSEA in multi-core environment.I used the 1ktools(https://github.com/cmap/l1ktools) tools to parse the .gctx file which stored gene profile data defined by Lincs(CMap) based on HDF5 file format.I used Matlab to parse the .gctx file、extract the gene profile sets and write to .txt file. C will read the file 、complete parallel GSEA and write out the result in a suitable file. 

Typically, in our work, we first implemented GSEA approach in efficient parallel strategy with MPI and OpenMP.In this part, on the one hand, we reduced the computational overhead of standard procedure to calculate the Enrichment Score by pre-sorting, indexing and removing the prefix sum. On the other hand, we will take a global permutation method to wipe off the redundant overhead of estimation of significance level step and multiple hypothesis testing step. 

Second, we expanded GSEA’s application to quickly compare two gene profile sets to get an Enrichment Score matrix of every gene profile pairs. In this part, in addition to using the previous optimization strategies, our implementation also allows to generate a second level of parallelization by creating several threads per MPI process.

Third, we clustered the gene profile based on the Enrichment Score matrix which we can get by the second part. In this part, Enrichment Score is served as the metric to measure the similarity between two gene profiles. We implemented a general clustering algorithm like K-Mediods which is an improved version of K-Means. The algorithm can quickly converge and then output the corresponding results.

I will update the tools as they become available.

### Matlab Tools: matlab_for_parse/

#### Requirements:

1. Matlab R2009a and above

#### Setting the MATLAB path:
Run `cd paraGSEA/matlab_for_parse` or enter the `pathtool` command, click "Add with Subfolders...", and select the directory `paraGSEA/matlab_for_parse`.

#### using by command line:
After you set the MATLAB path, you should enter "matlab" in shell to start matlab environment.
Then, you can set the input and output file path in runPreGSEA.m and enter "runPreGSEA" to parse the data.
Or, you can the using the matlab scripts below.

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

#### Description of Matlab Outputs

| File | Description | Format|
| ---- | ----------- |--------|
| data_for_test.txt | ranked profile file | first line : profile_number	profile_Length
											next profile_number lines : a ranked profile file included profile_Length elements |
| data_for_test_cid.txt | cid file | each line : a cid string corresponding to the last outputfile last profile_number lines |


### C Tools: src/

#### Requirements:

1. MPI and gcc compiler supports OpenMP

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

#### Note:
 * runParalESLinux.sh annotate a list of execution case of C tools. Removing the annotation, you can using it easily.
 * the details of parameter list of each C tools can be seen in runParalESLinux.sh or the TUTORIAL.


## The LINCS Dataset

The CMAP Cloud API offers programmatic access to annotations and perturbational signatures in [the LINCS L1000 dataset](http://lincscloud.org/) via a collection of HTTP-based RESTful web services. You can get the .gctx file stored gene profile data by the Website.