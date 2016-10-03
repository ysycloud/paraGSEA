<a name="runGenReference.doc"></a>
# runGenReference #

This shell script starts matlab environment and parses original data using the tools
provided by 1ktools project. some reference data files will be wrote out for the next GSEA analysis.
The detail is described in `tutorial`.

There are several parameters should be setted in command line example.

| parameter name | Parameter direction | Parameter function |
| -------------- | ------------------- | ------------------ |
| datasource | input | the file is based on HDF5 format with a gctx suffix stored gene profile data defined by Lincs(CMap) |
| gene_symbol_rhd | input | gene symbol field name of Lincs defined in this datasource |
| sample_conditions_chd | input | sample conditions field names of Lincs defined in this datasource |

An example of how to set these parameters is given in this script, and a 
.gctx file is provided in `data/` directory. you can easily change them 
in your need.

Three files will generate as reference data(`Gene_List.txt`, `Samples_Condition.txt`, `Samples_RowByteOffset.txt`) in `data/Reference` directory


