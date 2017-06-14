<a name="getReferences.doc"></a>
# getReferences #

This tool is used when the `rhd` and `chd` structs have been splited from some new datasets (mainly `.gctx`) of LINCS to be separate text files. 
As for this kind of datasets, we cannot make sure the value of each field just through the `gctx` files and then we cannot extract certain profiles from these datasets.
Obviously, we also cannot get reference data from these new datasets. Therefore, `genReferenceforNewDataSet.m` is replaced by this tool. Some reference data files will be wrote out for the next GSEA analysis.
The detail is described in `tutorial`.

There are several parameters should be setted in command line example.

| parameter name | Parameter function |
| -------------- | ------------------ |
| -1 --input1 | gene_info txt file |
| -2 --input2 | inst_info txt file |
| -r --reference | a directory used for outputing referenced files about genesymbols and cids |

An example of how to set these parameters and run this tool is given in `runGetReferencesLinux.sh` script. the examples of `gene_info` and `inst_info` txt file is provided in `data/GSE92742_INFO.tar.gz`. you can easily change them in your need.

Three files will generate as reference data(`Gene_List.txt`, `Samples_Condition.txt`, `Samples_RowByteOffset.txt`) in directory you set in `-r`. the example results are also stored in `data/GSE92742_INFO/Reference` after you decompress `data/GSE92742_INFO.tar.gz`.


