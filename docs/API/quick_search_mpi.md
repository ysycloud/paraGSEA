<a name="quick_search_mpi.doc"></a>
# quick_search_mpi #

This tool implemented original GSEA approach with GeneSet and Profiles
which will complete parallel GSEA in a multiprocesses way by MPI and 
show the topN results.

There are several parameters should be setted in command line example.

| parameter name | Parameter function |
| -------------- | -------------------|
| -n process_num | Total number of processes |
| -ppn pernum |the number of processes in each node |
| -hostfile hostfile | list the IP or Hostname of nodes |
| filename | a parsed profiles's file from pretreatment stage |
| topN | The first and last N GSEA records ordered by ES |

A sample `Shell script` file is given below that makes use of `quick_search_mpi`.

```shell
mpirun -n 2 -ppn 2 -hostfile example/hostfile ./bin/quick_search_mpi "data/data_for_test.txt" 15
```

the example hostfile just has one record as below:
```shell
localhost
```

However, there is a thing we must notice:
you'd better keep the formula `process_num = pernum * number of IP(Hostname) list in hostfile` 
correct.

It may will produce the following output:
```shell
Profile Set is Loading...!
profilenum:272	 genelen:978
loading IO and prework time by MPI: 0.0244 s
input the GeneSet( a integer[1-genelen] string split by space ):
```

Then you can input the GeneSet which is a integer string split by space 
between 1 to genlen. each integer represents a gene recorded in 
`data/data_for_test_rid.txt`

Therefore, you can input the GeneSet as below:
```shell
100 99 88 87 86 85 84 83 82 81
```

It will produce the following output:
```shell
printf the high level of TopN GSEA result:
NO.1 -> cid:204  ES:0.528926  NES:2.080132  pv:0.0016774794
NO.2 -> cid:66  ES:0.518595  NES:1.985430  pv:0.0016865702
NO.3 -> cid:37  ES:0.456818  NES:1.780453  pv:0.0059239669
NO.4 -> cid:78  ES:0.442355  NES:1.703747  pv:0.0098807840
NO.5 -> cid:97  ES:0.436364  NES:1.715039  pv:0.0072045455
NO.6 -> cid:255  ES:0.409711  NES:1.621115  pv:0.0101367769
NO.7 -> cid:65  ES:0.397107  NES:1.518632  pv:0.0160305767
NO.8 -> cid:206  ES:0.396488  NES:1.561015  pv:0.0147049580
NO.9 -> cid:58  ES:0.393595  NES:1.529173  pv:0.0147747936
NO.10 -> cid:4  ES:0.391942  NES:1.516823  pv:0.0148828506
NO.11 -> cid:184  ES:0.383884  NES:1.470641  pv:0.0172006207
NO.12 -> cid:183  ES:0.381818  NES:1.481119  pv:0.0159822311
NO.13 -> cid:71  ES:0.381198  NES:1.454726  pv:0.0210621910
NO.14 -> cid:185  ES:0.379752  NES:1.467693  pv:0.0200134277
NO.15 -> cid:123  ES:0.373967  NES:1.458159  pv:0.0213320255
printf the low level of TopN GSEA result:
NO.1 -> cid:137  ES:-0.613636  NES:-2.279210  pv:0.0000000000
NO.2 -> cid:91  ES:-0.602273  NES:-2.284496  pv:-0.0006931818
NO.3 -> cid:90  ES:-0.465083  NES:-1.796902  pv:-0.0071681814
NO.4 -> cid:77  ES:-0.446281  NES:-1.704311  pv:-0.0097188015
NO.5 -> cid:237  ES:-0.427066  NES:-1.622913  pv:-0.0089128094
NO.6 -> cid:268  ES:-0.423760  NES:-1.613431  pv:-0.0100869846
NO.7 -> cid:17  ES:-0.418388  NES:-1.555209  pv:-0.0177227268
NO.8 -> cid:124  ES:-0.377893  NES:-1.405311  pv:-0.0274190063
NO.9 -> cid:52  ES:-0.375826  NES:-1.398214  pv:-0.0236545448
NO.10 -> cid:60  ES:-0.374380  NES:-1.423725  pv:-0.0251080551
NO.11 -> cid:247  ES:-0.370041  NES:-1.397852  pv:-0.0228958721
NO.12 -> cid:89  ES:-0.367562  NES:-1.428832  pv:-0.0219681835
NO.13 -> cid:63  ES:-0.365289  NES:-1.379106  pv:-0.0240452442
NO.14 -> cid:135  ES:-0.359504  NES:-1.388710  pv:-0.0260508270
NO.15 -> cid:270  ES:-0.356818  NES:-1.336288  pv:-0.0287752094
finish GSEA time by MPI: 0.1705 s
input the GeneSet(split by space):
```

Then you can input the GeneSet again to continue to complete GSEA approach,
or, you can input `exit` to stop this tool as below:
```shell
input the GeneSet(split by space):
exit
```