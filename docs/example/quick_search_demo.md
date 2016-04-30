<a name="quick_search_demo.doc"></a>
# quick_search_demo #

This shell script uses `runPreGSEAbyMatlab.sh` or `runparaPreGSEAbyMatlab.sh` to parse original profile 
data first, then, three C Tools, which you can chose, of Quick_Search implemented original GSEA approach 
with GeneSet and Profiles and show the topN results will be used. 

There is one parameter you should notice to input when execute the script.

| parameter name | Parameter function |
| -------------- | ------------------ |
| quick_search_way | which Tool do you chose to complete the GSEA analysis(0_serial,1_openmp,2_mpi) |

which tools to be chosen is decided by your need. However, there are still some parameters should be setted 
in command line example when you use these tools. An example of how to set these parameters is given in this script. 
you can easily change them in your need.

The detailed usage of these tools is described in `docs\Tools\quick_search_*.md`.
