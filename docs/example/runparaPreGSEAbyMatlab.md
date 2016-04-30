<a name="runparaPreGSEAbyMatlab.doc"></a>
# runparaPreGSEAbyMatlab #

This shell script starts matlab environment and parses original data in a parallel way using the tools
provided by 1ktools project. some files will be wrote out for the next GSEA analysis.

There are several parameters should be setted in command line example.

| parameter name | Parameter direction | Parameter function |
| -------------- | ------------------- | ------------------ |
| file_input | input | the file is based on HDF5 format with a gctx suffix stored gene profile data defined by Lincs(CMap)|
| file_name | output | ranked profile file in a index format |
| file_name_cid | output | cid( profile identification ) file, each line is a cid( profile identification ) string  |
| file_name_rid | output | rid( gene identification ) file, each line is a rid( gene identification ) string corresponding to the rid attribute of input file and its index corresponding to the ranked profile file |
| cores | input | the number of cores you want to use in the parse work |

An example of how to set these parameters is given in this script, and a 
.gctx file is provided in `data/` directory. you can easily change them 
in your need.

However, there are some things we must notice:

1. you are supposed to make sure that you have a multicores system first.

2. the number of cores must be smaller than the actual core number in your system. 

3. after the parse work, you shoul merge every parts of output file into a whole file like the shell script shown.
```shell
cat ../data/data_for_test.txt_* >> ../data/data_for_test.txt
cat ../data/data_for_test_cid.txt_* >> ../data/data_for_test_cid.txt
rm -f ../data/data_for_test.txt_* ../data/data_for_test_cid.txt_*
```
therefore, in the merge work, the file name must be same with the `file_name` and `file_name_cid` when you started 
and executed the Matlab script.
```shell
# execute matlab script to parallel parse the data  
matlab -nodesktop -nosplash -nojvm -r "file_input='../data/modzs_n272x978.gctx'; file_name='../data/data_for_test.txt'; file_name_cid='../data/data_for_test_cid.txt'; file_name_rid='../data/data_for_test_rid.txt'; cores = 2; paraPreGSEA; quit;"
```

