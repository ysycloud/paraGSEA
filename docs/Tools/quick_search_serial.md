<a name="quick_search_serial.doc"></a>
# quick_search_serial #

This tool implemented original GSEA approach with GeneSet and Profiles
which will complete GSEA and show the topN results in a serial way.

There are several parameters should be setted in command line example.

| parameter name | Parameter function |
| -------------- | -------------------|
| -i --input | a parsed profiles's ranked sequence number file from pretreatment stage |
| -n --topn | define the first and last N GSEA records ordered by ES |
| -s --sample | a text file include sample sequence numbers which are extracted from pretreatment stage |
| -r --reference | a directory include some reference data files we generate from pretreatment stage |


A sample `Shell script` file is given below that makes use of `quick_search_serial`.

```shell
./bin/quick_search_serial -i data/data_for_test.txt -n 8 -s data/data_for_test_cidnum.txt -r data/Reference
```

It may will produce the following output:
```shell
Profile Set is Loading...!
profilenum:272	 genelen:978
loading IO and prework time: 0.0237 s
which way do you want to input the GeneSet( 0 -> standard input , others -> file input ):
```

Then you can choose which way do you want to input the GeneSet. 
Choosing 0, then you must input a gene set directly. 

For example

Choosing Others, you must input a file path where there is a gene set. 
Second way may be more convenient such as the example shows.


Therefore, you can input the GeneSet as below:
```shell
2 4 21 53 33 15 99 32 6 100
```

It will produce the following output:
```shell
printf the high level of TopN GSEA result:
NO.1 -> cid:70  ES:0.591942  NES:2.289738  pv:0.0000000000
NO.2 -> cid:41  ES:0.530165  NES:2.083023  pv:0.0005632231
NO.3 -> cid:207  ES:0.494008  NES:1.940730  pv:0.0020605371
NO.4 -> cid:34  ES:0.490909  NES:1.868850  pv:0.0062698340
NO.5 -> cid:67  ES:0.465289  NES:1.752064  pv:0.0051654959
NO.6 -> cid:185  ES:0.458884  NES:1.805547  pv:0.0043198342
NO.7 -> cid:256  ES:0.455165  NES:1.744763  pv:0.0052508264
NO.8 -> cid:82  ES:0.453099  NES:1.795171  pv:0.0062462807
NO.9 -> cid:71  ES:0.453099  NES:1.762074  pv:0.0062785120
NO.10 -> cid:129  ES:0.451653  NES:1.772830  pv:0.0049888430
printf the low level of TopN GSEA result:
NO.1 -> cid:112  ES:-0.499174  NES:-1.921487  pv:-0.0010840909
NO.2 -> cid:268  ES:-0.485537  NES:-1.912314  pv:-0.0014842975
NO.3 -> cid:237  ES:-0.457025  NES:-1.752636  pv:-0.0067721071
NO.4 -> cid:119  ES:-0.447521  NES:-1.735601  pv:-0.0040384297
NO.5 -> cid:267  ES:-0.439669  NES:-1.715914  pv:-0.0048855371
NO.6 -> cid:246  ES:-0.433058  NES:-1.673301  pv:-0.0048192153
NO.7 -> cid:42  ES:-0.429752  NES:-1.685708  pv:-0.0031320248
NO.8 -> cid:98  ES:-0.429339  NES:-1.638468  pv:-0.0102811985
NO.9 -> cid:236  ES:-0.417149  NES:-1.589109  pv:-0.0126154938
NO.10 -> cid:54  ES:-0.416322  NES:-1.575946  pv:-0.0116692152
finish GSEA time: 0.0864 s
input the GeneSet(split by space):
```

Then you can input the GeneSet again to continue to complete GSEA approach,
or, you can input `exit` to stop this tool as below:
```shell
input the GeneSet(split by space):
exit
```