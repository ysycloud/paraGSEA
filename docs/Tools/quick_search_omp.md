<a name="quick_search_omp.doc"></a>
# quick_search_omp #

This tool implemented original GSEA approach with GeneSet and Profiles which 
will complete parallel GSEA in a multithreads way by OpenMP and show the topN
results.

There are several parameters should be setted in command line example.

| parameter name | Parameter function |
| -------------- | -------------------|
| filename | a parsed profiles's file from pretreatment stage |
| thread_num | the number of threads |
| topN | The first and last N GSEA records ordered by ES |
| paraway | parallel ways( 0: split data in balance load way by ourselves; 1: using #pragma omp for ) |

A sample `Shell script` file is given below that makes use of `quick_search_omp`.

```shell
./bin/quick_search_omp "data/data_for_test.txt" 4 10 1
```

It may will produce the following output:
```shell
Profile Set is Loading...!
profilenum:272	 genelen:978
loading IO and prework time by openmp: 0.0267 s
input the GeneSet( a integer[1-genelen] string split by space ):
```

Then you can input the GeneSet which is a integer string split by space 
between 1 to genlen. each integer represents a gene recorded in 
`data/data_for_test_rid.txt`

Therefore, you can input the GeneSet as below:
```shell
1 2 3 4 5 6 7 8 9
```

It will produce the following output:
```shell
printf the high level of TopN GSEA result:
NO.1 -> cid:53  ES:0.543516  NES:1.961615  pv:0.0011159270
NO.2 -> cid:70  ES:0.513932  NES:1.900128  pv:0.0021809425
NO.3 -> cid:269  ES:0.495356  NES:1.832136  pv:0.0020577917
NO.4 -> cid:43  ES:0.488132  NES:1.811714  pv:0.0047254906
NO.5 -> cid:196  ES:0.470932  NES:1.764590  pv:0.0015411077
NO.6 -> cid:209  ES:0.437909  NES:1.578733  pv:0.0123808050
NO.7 -> cid:240  ES:0.427245  NES:1.563762  pv:0.0123505325
NO.8 -> cid:191  ES:0.422773  NES:1.580992  pv:0.0153061571
NO.9 -> cid:184  ES:0.421053  NES:1.553085  pv:0.0157031298
NO.10 -> cid:20  ES:0.404541  NES:1.470701  pv:0.0152803593
printf the low level of TopN GSEA result:
NO.1 -> cid:143  ES:-0.503956  NES:-1.838501  pv:-0.0037485380
NO.2 -> cid:258  ES:-0.493636  NES:-1.793665  pv:-0.0065067077
NO.3 -> cid:236  ES:-0.492948  NES:-1.773078  pv:-0.0038407295
NO.4 -> cid:237  ES:-0.436189  NES:-1.571809  pv:-0.0102404547
NO.5 -> cid:44  ES:-0.434469  NES:-1.571682  pv:-0.0139305115
NO.6 -> cid:226  ES:-0.433093  NES:-1.582485  pv:-0.0113460636
NO.7 -> cid:112  ES:-0.428621  NES:-1.533743  pv:-0.0134251823
NO.8 -> cid:94  ES:-0.428277  NES:-1.539260  pv:-0.0158121767
NO.9 -> cid:250  ES:-0.416581  NES:-1.503345  pv:-0.0161757813
NO.10 -> cid:254  ES:-0.413485  NES:-1.509609  pv:-0.0135163393
finish GSEA time by openmp: 0.1124 s
input the GeneSet(split by space):
```

Then you can input the GeneSet again to continue to complete GSEA approach,
or, you can input `exit` to stop this tool as below:
```shell
input the GeneSet(split by space):
exit
```