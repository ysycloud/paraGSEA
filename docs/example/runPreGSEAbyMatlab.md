<a name="runPreGSEAbyMatlab.doc"></a>
# runPreGSEAbyMatlab #

This shell script starts matlab environment and parses original data using the tools
provided by 1ktools project. some files will be wrote out for the next GSEA analysis.

There are several parameters should be setted in command line example.

| parameter name | Parameter direction | Parameter function |
| -------------- | ------------------- | ------------------ |
| file_input | input | the file is based on HDF5 format with a gctx suffix stored gene profile data defined by Lincs(CMap)|
| file_name | output | ranked profile file in a index format |
| file_name_cid | output | cid( profile identification ) file, each line is a cid( profile identification ) string  |
| file_name_rid | output | rid( gene identification ) file, each line is a rid( gene identification ) string corresponding to the rid attribute of input file and its index corresponding to the ranked profile file |

A example of how to set these parameters is given in this script, and a 
.gctx file is provided in `data/` directory. you can easily change them 
in your need.


