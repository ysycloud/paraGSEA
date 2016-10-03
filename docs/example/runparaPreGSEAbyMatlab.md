<a name="runparaPreGSEAbyMatlab.doc"></a>
# runparaPreGSEAbyMatlab #

This shell script starts matlab environment and parses original data in a parallel way using the tools
provided by 1ktools project. some files will be wrote out for the next GSEA analysis.
The detail is described in `tutorial`.

There are several parameters should be setted in command line example.

| parameter name | Parameter direction | Parameter function |
| -------------- | ------------------- | ------------------ |
| file_input | input | the file is based on HDF5 format with a gctx suffix stored gene profile data defined by Lincs(CMap)|
| file_name | output | ranked profiles file in a sequence number format |
| file_name_cidnum | output | profile sequence number file corresponding to these profiles we extract in pretreatment stage |
| cores | input | the number of cores you want to use in the parse work |
| sample_conditions_chd | input | sample conditions field names of Lincs defined in this datasource |
| cell_id_set | input | cell lines set where the profiles should be get from |
| pert_set | input | perturbations that should be used in experiments to get the profiles |
| pert_type_set | input | perturbation types that should be used in experiments to get the profiles |
| duration | input | duration to carry out the experiments |
| concentration | input |  concentration is supposed to be kept during the experiments |

An example of how to set these parameters is given in this script, and a 
.gctx file is provided in `data/` directory. you can easily change them 
in your need.

However, there are some things we must notice:

1. you are supposed to make sure that you have a multicores system first.

2. the number of cores must be smaller than the actual core number in your system.

3. the five parameters is the filter conditions to extract the profiles we need. 
If there is one condition we not set, this condition will no longer be taken into account to extract profiles. 

4. after the parse work, you shoul merge every parts of output file into a whole file like the shell script shown.
```shell
cat ../data/data_for_test.txt_* >> ../data/data_for_test.txt
cat ../data/data_for_test_cid.txt_* >> ../data/data_for_test_cid.txt
rm -f ../data/data_for_test.txt_* ../data/data_for_test_cid.txt_*
```
therefore, in the merge work, the file name must be same with the `file_name` and `file_name_cid` when you started 
and executed the Matlab script.

