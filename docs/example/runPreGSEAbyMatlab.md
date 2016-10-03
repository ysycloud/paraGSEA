<a name="runPreGSEAbyMatlab.doc"></a>
# runPreGSEAbyMatlab #

This shell script starts matlab environment and parses original data using the tools
provided by 1ktools project. some files will be wrote out for the next GSEA analysis.
The detail is described in `tutorial`.

There are several parameters should be setted in command line example.

| parameter name | Parameter direction | Parameter function |
| -------------- | ------------------- | ------------------ |
| file_input | input | the file is based on HDF5 format with a gctx(gct) suffix stored gene profile data defined by Lincs(CMap)|
| file_name | output | ranked profiles file in a sequence number format |
| file_name_cidnum | output | profile sequence number file corresponding to these profiles we extract in pretreatment stage |
| sample_conditions_chd | input | sample conditions field names of Lincs defined in this datasource |
| cell_id_set | input | cell lines set where the profiles should be get from |
| pert_set | input | perturbations that should be used in experiments to get the profiles |
| pert_type_set | input | perturbation types that should be used in experiments to get the profiles |
| duration | input | duration to carry out the experiments |
| concentration | input |  concentration is supposed to be kept during the experiments |

An example of how to set these parameters is given in this script, and a 
.gctx file is provided in `data/` directory. you can easily change them 
in your need.

Note that the five parameters is the filter conditions to extract the profiles we need.
If there is one condition we not set, this condition will no longer be taken into account to extract profiles. 

Also, after the parse work, you shoul merge some intermediate files into a whole file like the shell script shown.
```shell
cat ../data/tmp >> ../data/data_for_test.txt
rm -f ../data/tmp
```

