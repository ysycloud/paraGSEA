<a name="cluster_demo.doc"></a>
# cluster_demo #

This shell script uses `runPreGSEAbyMatlab.sh` or `runparaPreGSEAbyMatlab.sh` to parse original profile 
data first, then, three C Tools of ES_Matrix which compare two gene profile sets to get an Enrichment Score matrix of every gene profile pairs are supposed to be chosen to perform and write out an ES Matrix in a distributed way. Finally, two C Tools of cluster using ES Matrix as an input are also supposed to be chosen to perform and write out clustering results in the root node.

There are two parameters you should notice to input when execute the script.

| parameter name | Parameter function |
| -------------- | ------------------ |
| matrix_way | which Tool do you chose to get ES Matrix(0_nocom,1_p2p,2_cocom[when distribute the input file]) |
| cluster_way | which Tool do you chose to complete the Cluster algorithm(0_KMediods,1_KMediods++) |

which tools to be chosen is decided by your need. However, there are still some parameters should be setted 
in command line example when you use these tools. An example of how to set these parameters is given in this script. 
you can easily change them in your need.

The detailed usage of these tools is described in `docs\Tools\ES_Matrix_ompi_*.md` and `docs\Tools\Cluster_*_ompi.md`.
Also, you can see the detail usage in `tutorial`.

However, there are still some things you must notice when you set these parameters:

1. you'd better keep the formula `process_num = pernum * number of IP(Hostname) list in hostfile` correct.

2. you must note that input matrix should include the same identity of rows and columns, which means the program that calculates ES Matrix is supposed to use same two file as its input. Only in this way can we get the similarity of each profile pair.

3. you must also note that the MPI Settings and hostfile should not be changed compared to the program that calculates ES Matrix. Because the ES_Matrix is stored in distributed way, if you change these settings, each process can not find the right ES matrix blocks.

The sample script have guaranteed these condition above, you are not supposed to change these rules when you change these parameters.



