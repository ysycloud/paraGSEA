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
| -i --input | a parsed profiles's ranked sequence number file from pretreatment stage |
| -n --topn | define the first and last N GSEA records ordered by ES |
| -s --sample | a text file include sample sequence numbers which are extracted from pretreatment stage |
| -r --reference | a directory include some reference data files we generate from pretreatment stage |

A sample `Shell script` file is given below that makes use of `quick_search_mpi`.

```shell
mpirun -n 2 -ppn 2 -hostfile example/hostfile quick_search_mpi -i data/data_for_test.txt -n 5 -s data/data_for_test_cidnum.txt -r data/Reference
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
which way do you want to input the GeneSet( 0 -> standard input , others -> file input ):
```

Then you can choose which way do you want to input the GeneSet. 
Choosing 0, then you must input a gene set directly. 

For example:
```shell
input the GeneSet until 'exit'( a string of each Gene Symbol split by space ):
CCNH HMGA2 IGFBP3 RB1 PARP1 CDK6
```

Or choosing Others, you must input a file path where there is a gene set. 
Most of times, Second way may be more convenient.

For example:
```shell
input the path of file that has GeneSet until 'exit'(each line has a Gene Symbol/name):
data/GeneSet.txt
```

It will produce the following output:
```shell

printf the high level of TopN GSEA result:

NO.1 -> SampleConditions: cid:CPC006_SKLU1_6H:BRD-K56343971-001-02-3:10;    cell_line:     SKLU1;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration  ES:0.315086  NES:2.540525  pv:0.0000000000

NO.2 -> SampleConditions: cid:LJP001_MCF10A_24H:BRD-K56343971-001-04-9:0.08;    cell_line:    MCF10A;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    dura  ES:0.225345  NES:1.881119  pv:0.0010856895

NO.3 -> SampleConditions: cid:NMH001_NEU.KCL_6H.4H:BRD-K69726342-001-02-6:10;    cell_line:   NEU.KCL;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    dur  ES:0.223448  NES:1.843819  pv:0.0014393966

NO.4 -> SampleConditions: cid:LJP001_BT20_24H:BRD-K56343971-001-04-9:0.4;    cell_line:      BT20;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duratio  ES:0.221078  NES:1.843653  pv:0.0016556464

NO.5 -> SampleConditions: cid:LJP001_MCF7_24H:BRD-K56343971-001-04-9:2;    cell_line:      MCF7;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:  ES:0.216035  NES:1.790939  pv:0.0023914223

printf the low level of TopN GSEA result:

NO.1 -> SampleConditions: cid:LJP001_HS578T_6H:BRD-K56343971-001-04-9:0.4;    cell_line:    HS578T;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    durati  ES:-0.229569  NES:-1.821155  pv:-0.0023933623 

NO.2 -> SampleConditions: cid:CPC006_HEC108_6H:BRD-U88459701-000-01-8:10;    cell_line:    HEC108;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duratio  ES:-0.219655  NES:-1.762991  pv:-0.0029483190

NO.3 -> SampleConditions: cid:CPC006_SNUC5_6H:BRD-K56343971-001-02-3:10;    cell_line:     SNUC5;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration  ES:-0.202629  NES:-1.648118  pv:-0.0058889213

NO.4 -> SampleConditions: cid:CPC004_A375_6H:BRD-A51714012-001-03-1:10;    cell_line:      A375;    perturbation:    venlafaxine;    perturbation type:    trt_cp;    duration:  ES:-0.194224  NES:-1.559712  pv:-0.0068219395

NO.5 -> SampleConditions: cid:CPC006_SNGM_6H:BRD-U88459701-000-01-8:10;    cell_line:      SNGM;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:  ES:-0.192112  NES:-1.541771  pv:-0.0081038799

finish GSEA time by MPI: 0.1705 s
```

Then you can input the GeneSet again to continue to complete GSEA approach,
or, you can input `exit` to stop this tool as below:
```shell
input the path of file that has GeneSet until 'exit'(each line has a Gene Symbol/name):
exit
```