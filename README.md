# ParaGSEA
Gene set Enrichment analysis(GSEA) Tools for Lincs data in a parallel manner

## Tools

A brief description of the tools is given below. The tools implement the function of parallel compute GSEA in multi-core environment.I used the 1ktools(https://github.com/cmap/l1ktools) tools to parse the .gctx file which stored gene profile data defined by Lincs(CMap) based on HDF5 file format.I used Matlab to parse the .gctx file、extract the gene profile sets and write to .txt file. C will read the file 、complete parallel GSEA and write out the result in a suitable file. 
Typically, in our work, we first implemented GSEA approach in efficient parallel strategy with MPI and OpenMP.In this part, on the one hand, we reduced the computational overhead of standard procedure to calculate the Enrichment Score by pre-sorting, indexing and removing the prefix sum. On the other hand, we will take a global permutation method to wipe off the redundant overhead of estimation of significance level step and multiple hypothesis testing step. 
Second, we expanded GSEA’s application to quickly compare two gene profile sets to get an Enrichment Score matrix of every gene profile pairs. In this part, in addition to using the previous optimization strategies, our implementation also allows to generate a second level of parallelization by creating several threads per MPI process.
Third, we clustered the gene profile based on the Enrichment Score matrix which we can get by the second part. In this part, Enrichment Score is served as the metric to measure the similarity between two gene profiles. We implemented a general clustering algorithm like K-Means. The algorithm can quickly converge and then output the corresponding results.
I will update the tools as they become available.

### Matlab Tools: matlab/

#### Requirements:

1. Matlab R2009a and above

#### Setting the MATLAB path:
Enter the "pathtool" command, click "Add with Subfolders...", and select the directory ParaGSEA/matlab_for_parse.

#### using by command line:
After you set the MATLAB path, you should enter "matlab" in shell to start matlab environment.
Then, you can set the input and output file path in runPreGSEA.m and enter "runPreGSEA" to parse the data.

#### Tools:
* [**runPreGSEA.m**] : Setting Parameters and execute the PreGSEA.
* [**PreGSEA.m**] : extract the gene profile sets、finish pre-sorting and write to .txt file.


### C Tools: c/

#### Requirements:

1. MPI and gcc compiler supports OpenMP

#### Tools:

* [**quick_search_serial.c**] read the .txt file 、complete GSEA and show the topN results in a serial way.
* [**quick_search_omp.c**] read the .txt file 、complete parallel GSEA by OpenMP and show the topN results.
* [**quick_search_mpi.c**] read the .txt file 、complete parallel GSEA by MPI and show the topN results.
* [**ES_Matrix_ompi_nocom.c**] read the .txt file 、complete parallel computing ES_Matrix and write out the result by MPI/OpenMP with no communication in distributing second file of input.
* [**ES_Matrix_ompi_p2p.c**] read the .txt file 、complete parallel computing ES_Matrix and write out the result by MPI/OpenMP with p2p communication in distributing second file of input.
* [**ES_Matrix_ompi_cocom.c**] read the .txt file 、complete parallel computing ES_Matrix and write out the result by MPI/OpenMP with collective communication in distributing second file of input.
* [**Cluster_KMeans_ompi.c**] read the ES_Matrix file 、complete a general clustering algorithm like K-Means by MPI/OpenMP.

#### Demo:
* [**compileParalES.sh**]: compile the C source code to Executable files.
* [**runParalESLinux.sh**]: run Executable files in Linux.

#### Note:
 * compileParalES.sh will carry out a makefile file which is finished in current directory.
 * runParalESLinux.sh annotate a list of execution case of C tools. Removing the annotation, you can using it easily.
 * the parameter list of each C tools can be seen in runParalESLinux.sh


## The LINCS Dataset

The CMAP Cloud API offers programmatic access to annotations and perturbational signatures in [the LINCS L1000 dataset](http://lincscloud.org/) via a collection of HTTP-based RESTful web services. You can get the .gctx file stored gene profile data by the Website.