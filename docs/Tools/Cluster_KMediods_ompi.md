<a name="Cluster_KMediods_ompi.doc"></a>
# Cluster_KMediods_ompi #

This tool read the ES_Matrix file in a distributed way and 
cluster the gene profiles based on the Enrichment Score matrix 
which can be seemed as the similarity Matrix of gene profiles
by MPI/OpenMP.

There are several parameters should be setted in command line example.

| parameter name | Parameter function |
| -------------- | -------------------|
| -n process_num | Total number of processes |
| -ppn pernum |the number of processes in each node |
| -hostfile hostfile | list the IP or Hostname of nodes |
| -t --thread | define the maximum number of parallel threads |
| -c --cluster | define the number of clusters we want to get |
| -w --write | decide whether output the results |
| -i --input | distributed ES_Matrix file we get from stage 2(Compare Profiles) |
| -o --output | define the output cluster result file of every profiles in root node |
| -s --sample | a text file include sample sequence numbers which are extracted from pretreatment stage |
| -r --reference | a directory include some reference data files we generate from pretreatment stage |

A sample `Shell script` file is given below that makes use of `Cluster_KMediods_ompi`.

```shell
mpirun -n 2 -ppn 2 -hostfile example/hostfile Cluster_KMediods_ompi -t 4 -c 8 -w 1 -i data/ES_Matrix_test -o data/Cluster_result_test.txt  -s data/data_for_test_cidnum.txt -r data/Reference
```

the example hostfile just has one record as below:
```shell
localhost
```

However, there are some things we must notice:

1. you'd better keep the formula `process_num = pernum * number of IP(Hostname) list in hostfile` correct.

2. When we want to execute the Cluster operator, we must note that input matrix should include the same identity of rows and columns, which means the program that calculates ES Matrix is supposed to use same two file as its input. Only in this way can we get the similarity of each profile pair.

3. When we want to execute the Cluster operator, we must also note that the MPI Settings and hostfile should not be changed compared to the program that calculates ES Matrix. Because the ES_Matrix is stored in distributed way, if you change these settings, each process can not find the right ES matrix blocks.

4. If you set the number of clusters too big, clustering algorithm may not converge quickly.


It may produce the following output:
```shell
Matrix is Loading...!
loading IO and prework time : 0.0214 s
Paral KMediods compute the Cluster Centers is Starting...!
Init cluster centers is:
127 216 145 245 71 29 265 35 
1th iteraction cluster_center_new is:
174 30 136 60 86 58 265 38 
2th iteraction cluster_center_new is:
18 198 68 120 0 136 212 265 
3th iteraction cluster_center_new is:
63 30 42 248 136 88 140 265 
4th iteraction cluster_center_new is:
18 30 63 122 136 208 240 265 
5th iteraction cluster_center_new is:
56 60 63 220 136 46 240 265 
6th iteraction cluster_center_new is:
20 56 0 200 136 248 42 265 
7th iteraction cluster_center_new is:
63 128 30 112 220 200 248 265 
8th iteraction cluster_center_new is:
38 56 66 174 0 192 248 265 
9th iteraction cluster_center_new is:
63 76 112 136 174 220 32 265 
10th iteraction cluster_center_new is:
231 200 76 256 136 174 248 265 
11th iteraction cluster_center_new is:
30 192 174 0 56 104 256 265 
12th iteraction cluster_center_new is:
63 30 112 208 212 248 22 265 
13th iteraction cluster_center_new is:
66 60 63 136 46 220 164 265 
14th iteraction cluster_center_new is:
86 0 128 132 174 248 140 265 
15th iteraction cluster_center_new is:
63 130 220 136 152 212 164 265 
16th iteraction cluster_center_new is:
63 156 136 152 248 222 184 264 
17th iteraction cluster_center_new is:
63 136 144 176 184 30 248 264 
18th iteraction cluster_center_new is:
18 63 136 208 164 184 84 264 
19th iteraction cluster_center_new is:
56 128 0 136 248 184 260 264 
20th iteraction cluster_center_new is:
63 56 220 136 184 104 260 264 
21th iteraction cluster_center_new is:
218 63 248 136 184 148 262 264 
22th iteraction cluster_center_new is:
63 136 208 184 30 248 262 264 
23th iteraction cluster_center_new is:
18 63 220 184 266 248 262 264 
24th iteraction cluster_center_new is:
56 63 174 220 38 262 264 102 
25th iteraction cluster_center_new is:
136 56 63 270 174 248 262 264 
26th iteraction cluster_center_new is:
140 63 220 174 248 262 264 270 
27th iteraction cluster_center_new is:
63 140 174 192 248 262 264 270 
28th iteraction cluster_center_new is:
63 140 174 192 248 262 264 270 
Paral KMediods	compute the Cluster Centers Spent: 0.2829 s
```

If the `-w` parameter is 0, it will output:
```shell
Just run for test, no results output
```
before
```shell
Paral KMediods compute the Cluster Centers Spent: 0.2829 s
```

There is also no more need of you to input anything in command line. However,
the cluster results will be written to file `data/Cluster_result_test.txt` 
in root process if the ‘-w’ parameter is 0.

For Example:

	cluster 1 :
	cid:CPC006_A549_6H:BRD-U88459701-000-01-8:10;    cell_line:      A549;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	cid:CPC020_A375_6H:BRD-A82307304-001-01-8:10;    cell_line:      A375;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	......
	cluster 2 :
	cid:CPC020_HT29_6H:BRD-A82307304-001-01-8:10;    cell_line:      HT29;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	cid:CPC006_A375_6H:BRD-U88459701-000-01-8:10;    cell_line:      A375;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	......
	cluster 3 :
	cid:CPC006_A549_24H:BRD-U88459701-000-01-8:10;    cell_line:      A549;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:      24 h;    concentration:     10 uM
	cid:CPC006_AGS_6H:BRD-U88459701-000-01-8:10;    cell_line:       AGS;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	......
	cluster 4 :
	cid:CPC006_HA1E_6H:BRD-U88459701-000-01-8:10;    cell_line:      HA1E;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	cid:CPC006_HCC15_6H:BRD-U88459701-000-01-8:10;    cell_line:     HCC15;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	......
	cluster 5 :
	cid:CPC006_TYKNU_6H:BRD-K56343971-001-02-3:10;    cell_line:     TYKNU;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	cid:CPC006_VCAP_6H:BRD-K56343971-001-02-3:10;    cell_line:      VCAP;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	......
	cluster 6 :
	cid:CPC006_HEC108_6H:BRD-U88459701-000-01-8:10;    cell_line:    HEC108;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	cid:CPC006_SKLU1_6H:BRD-U88459701-000-01-8:10;    cell_line:     SKLU1;    perturbation:   atorvastatin;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	......
	cluster 7 :
	cid:CPC004_HCC515_24H:BRD-A51714012-001-03-1:10;    cell_line:    HCC515;    perturbation:    venlafaxine;    perturbation type:    trt_cp;    duration:      24 h;    concentration:     10 uM
	cid:CPC011_A549_6H:BRD-A51714012-003-09-4:10;    cell_line:      A549;    perturbation:    venlafaxine;    perturbation type:    trt_cp;    duration:       6 h;    concentration:     10 uM
	......
	cluster 8 :
	cid:BRAF001_A375_24H:BRD-U73308409-000-01-9:0.625;    cell_line:      A375;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:      24 h;    concentration:    500 nM
	cid:BRAF001_A375_6H:BRD-U73308409-000-01-9:0.625;    cell_line:      A375;    perturbation:    vemurafenib;    perturbation type:    trt_cp;    duration:       6 h;    concentration:    500 nM
	......