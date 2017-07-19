<a name="quick_search_profile.doc"></a>
# quick_search_profile #

This tool implemented parallel Enrichment score calculation by MPI/OpenMP between Profiles
where users should input a profile file path and the tool will search the library and finally 
show the topN results.

There are several parameters should be setted in command line example.

| parameter name | Parameter function |
| -------------- | -------------------|
| -n process_num | Total number of processes |
| -ppn pernum |the number of processes in each node |
| -hostfile hostfile | list the IP or Hostname of nodes |
| -t --thread | define the maximum number of parallel threads |
| -l--siglen | define the length of Gene Expression Signature |
| -n --topn | define the first and last N GSEA records ordered by ES |
| -i --input | a parsed profiles's ranked sequence number file from pretreatment stage as profile library|
| -s --sample | a text file include sample sequence numbers which are extracted from pretreatment stage |
| -r --reference | a directory include some reference data files we generate from pretreatment stage |

A sample `Shell script` file is given below that makes use of `quick_search_profile`.

```shell
mpirun -n 2 -ppn 2 -hostfile example/hostfile bin/quick_search_profile -i data/data_for_test.txt -l 50 -t 4 -n 8 -s data/data_for_test_cidnum.txt 
```

the example hostfile just has one record as below:
```shell
localhost
```

However, there is a thing we must notice:
you'd better keep the formula `process_num = pernum * number of IP(Hostname) list in hostfile` 
correct.

It may produce the following output:
```shell
Profile Set is Loading...!
profilenum:272   genelen:978
Memory check......
Available Memory:      260017644 KB
Needed Memory:      525 KB
loading IO and prework time: 0.0046 s
which type of profile do you want to input ( 0 -> sorted gene symbol list, others -> gene symbols and their expression levels ):
```

If Available Memory is less than Needed Memory, it will produce the following output:
```shell
available memory is not enough!!! Please use more nodes!!!
```
and finish the program.

If Available Memory is enough, you can choose which type of profile do you want to input.
Choosing 0, then you must input a file path of sorted gene symbol list. 

For example:
```shell
input the path of file that has target profile until 'exit'(each line has a Gene Symbol/name):
data/Profile.txt
```

The content of data/Profile.txt is as follows:
```shell
PSME1
ATF1
RHEB
FOXO3
RHOA
IL1B
ASAH1
RALA
ARHGEF12
SOX2
SERPINE1
HLA-DMA
EGF
APP
......
```

It will produce the following output:
```shell
printf the high level of TopN similar profile result:

NO.1 -> SampleConditions: cid:CPC006_TYKNU_6H:BRD-K56343971-001-02-3:10;    cell_line:     TYKNU;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:0.251045

NO.2 -> SampleConditions: cid:CPC011_MCF7_6H:BRD-A51714012-003-09-4:10;    cell_line:      MCF7;    perturbation:    venlafaxine;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:0.226487

NO.3 -> SampleConditions: cid:CPC004_HCC515_6H:BRD-A51714012-001-03-1:10;    cell_line:    HCC515;    perturbation:    venlafaxine;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:0.198847

NO.4 -> SampleConditions: cid:CPC006_MDST8_6H:BRD-U88459701-000-01-8:10;    cell_line:     MDST8;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:0.196067

NO.5 -> SampleConditions: cid:LJP001_HS578T_24H:BRD-K56343971-001-04-9:10;    cell_line:    HS578T;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:      24 h;    concentration:     10 ?M
  ES:0.196045

NO.6 -> SampleConditions: cid:LJP001_BT20_24H:BRD-K56343971-001-04-9:0.4;    cell_line:      BT20;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:      24 h;    concentration:    500 nM
  ES:0.186950

NO.7 -> SampleConditions: cid:LJP001_MCF10A_24H:BRD-K56343971-001-04-9:0.08;    cell_line:    MCF10A;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:      24 h;    concentration:    100 nM
  ES:0.186045

NO.8 -> SampleConditions: cid:CVD001_PHH_6H:BRD-U88459701-000-01-8:2.5;    cell_line:       PHH;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:      3 ?M
  ES:0.184946

printf the low level of TopN similar profile result:

NO.1 -> SampleConditions: cid:CPC004_A375_6H:BRD-A51714012-001-03-1:10;    cell_line:      A375;    perturbation:    venlafaxine;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:-0.338901

NO.2 -> SampleConditions: cid:CPC006_CORL23_6H:BRD-K56343971-001-02-3:10;    cell_line:    CORL23;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:-0.257672

NO.3 -> SampleConditions: cid:CPC006_RMUGS_6H:BRD-U88459701-000-01-8:10;    cell_line:     RMUGS;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:-0.232640

NO.4 -> SampleConditions: cid:CPC006_A549_24H:BRD-U88459701-000-01-8:10;    cell_line:      A549;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:      24 h;    concentration:     10 ?M
  ES:-0.209623

NO.5 -> SampleConditions: cid:CPC006_SW620_6H:BRD-U88459701-000-01-8:10;    cell_line:     SW620;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:-0.187112

NO.6 -> SampleConditions: cid:BRAF001_HEK293T_24H:BRD-U73308409-000-01-9:2.5;    cell_line:   HEK293T;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:      24 h;    concentration:      3 ?M
  ES:-0.185194

NO.7 -> SampleConditions: cid:CPC014_PC3_24H:BRD-K56343971-001-01-5:10;    cell_line:       PC3;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:      24 h;    concentration:     10 ?M
  ES:-0.166207

NO.8 -> SampleConditions: cid:CPC006_HCC15_6H:BRD-U88459701-000-01-8:10;    cell_line:     HCC15;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:-0.164774
finish ES time: 0.0023 s
input the path of file that has target profile until 'exit'(each line has a Gene Symbol/name):
```

If this file do not have the whole gene symbol in the library, there will be some notices and you should input again:
```shell
This profile is too short, please make sure it has same 978 genes with L1000 library profiles!
input the path of file that has target profile until 'exit'(each line has a Gene Symbol/name):
```


If you choosing Others at first, you must input a file path of gene symbols and their expression levels (spilted by `\t`). 

If you still input `data/Profile` without expression levels, the expression levels will be regarded as `0` and get the same results.

And If you input the file with correct format, take `data/ProfilewithExpression.txt` as an example ,
```shell
PSME1	0.45252848
ATF1	-0.018624753
RHEB	0.32816789
FOXO3	1.4956889
RHOA	1.5585777
IL1B	-0.39951164
ASAH1	-0.28378099
RALA	-0.63659459
ARHGEF12	-0.11138941
SOX2	-0.62310797
SERPINE1	-0.38765121
HLA-DMA	1.695841
EGF	-0.68614239
APP	-0.10260005
NOS3	-0.29481113
CSNK1A1	-0.15827303
NFATC4	-1.1765727
TBP	0.84266657
BRCA1	0.86614233
.......
```

It will produce the following output:
```shell

printf the high level of TopN similar profile result:

NO.1 -> SampleConditions: cid:CPC020_A375_6H:BRD-A82307304-001-01-8:10;    cell_line:      A375;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:1.000000

NO.2 -> SampleConditions: cid:CPC006_H1299_6H:BRD-U88459701-000-01-8:10;    cell_line:     H1299;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:0.391735

NO.3 -> SampleConditions: cid:NMH001_NEU_6H:BRD-K69726342-001-02-6:10;    cell_line:       NEU;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:0.389127

NO.4 -> SampleConditions: cid:CPC015_ASC_24H:BRD-A51714012-001-03-1:10;    cell_line:       ASC;    perturbation:    venlafaxine;    perturbation type:    trt_cp;    duration:      24 h;    concentration:     10 ?M
  ES:0.345140

NO.5 -> SampleConditions: cid:CPC015_HT29_6H:BRD-A51714012-001-03-1:10;    cell_line:      HT29;    perturbation:    venlafaxine;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:0.328588

NO.6 -> SampleConditions: cid:CPC006_VCAP_24H:BRD-U88459701-000-01-8:10;    cell_line:      VCAP;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:      24 h;    concentration:     10 ?M
  ES:0.324472

NO.7 -> SampleConditions: cid:LJP001_MCF10A_6H:BRD-K56343971-001-04-9:0.08;    cell_line:    MCF10A;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:       6 h;    concentration:    100 nM
  ES:0.322726

NO.8 -> SampleConditions: cid:CPC006_OV7_6H:BRD-K56343971-001-02-3:10;    cell_line:       OV7;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:0.320981

printf the low level of TopN similar profile result:

NO.1 -> SampleConditions: cid:CPC011_A549_6H:BRD-A51714012-003-09-4:10;    cell_line:      A549;    perturbation:    venlafaxine;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:-0.358966

NO.2 -> SampleConditions: cid:CPC006_HCT116_6H:BRD-K56343971-001-02-3:10;    cell_line:    HCT116;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:-0.344235

NO.3 -> SampleConditions: cid:CPC004_HEPG2_6H:BRD-A51714012-001-03-1:10;    cell_line:     HEPG2;    perturbation:    venlafaxine;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:-0.328772

NO.4 -> SampleConditions: cid:CPC004_HCC515_24H:BRD-A51714012-001-03-1:10;    cell_line:    HCC515;    perturbation:    venlafaxine;    perturbation type:    trt_cp;    duration:      24 h;    concentration:     10 ?M
  ES:-0.322489

NO.5 -> SampleConditions: cid:CPC006_RMUGS_6H:BRD-U88459701-000-01-8:10;    cell_line:     RMUGS;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:-0.292920

NO.6 -> SampleConditions: cid:LJP001_MDAMB231_6H:BRD-K56343971-001-04-9:0.08;    cell_line:  MDAMB231;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:       6 h;    concentration:    100 nM
  ES:-0.286509

NO.7 -> SampleConditions: cid:CPC014_VCAP_6H:BRD-K56343971-001-01-5:10;    cell_line:      VCAP;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 ?M
  ES:-0.275593

NO.8 -> SampleConditions: cid:LJP001_HS578T_6H:BRD-K56343971-001-04-9:0.4;    cell_line:    HS578T;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:       6 h;    concentration:    500 nM
  ES:-0.275280
finish ES time: 0.0038 s
```

Then you can input the file path again to continue,
or, you can input `exit` to stop this tool as below:
```shell
input the path of file that has target profile until 'exit'(each line has a Gene Symbol/name):
exit
```