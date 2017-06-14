<a name="Cluster_KMediods++_ompi.doc"></a>
# Cluster_KMediods++_ompi #

This tool read the ES_Matrix file in a distributed way and 
cluster the gene profiles based on the Enrichment Score matrix 
which can be seemed as the similarity Matrix of gene profiles
by MPI/OpenMP. However, let the distance between initial cluster
centers as far as possible

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

A sample `Shell script` file is given below that makes use of `Cluster_KMediods++_ompi`.

```shell
mpirun -n 2 -ppn 2 -hostfile example/hostfile Cluster_KMediods++_ompi -t 4 -c 8 -w 1 -i data/ES_Matrix_test -o data/Cluster_result_test.txt  -s data/data_for_test_cidnum.txt -r data/Reference
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
loading IO and prework time : 0.0168 s
Paral KMediods++ compute the Cluster Centers is Starting...!
Init cluster centers is:
23 76 194 150 138 1 11 15 
1th iteraction cluster_center_new is:
114 76 252 164 220 231 106 30 
2th iteraction cluster_center_new is:
18 76 178 240 248 184 136 0 
3th iteraction cluster_center_new is:
63 30 270 136 220 184 92 248 
4th iteraction cluster_center_new is:
18 126 28 136 184 208 248 268 
5th iteraction cluster_center_new is:
6 56 220 136 170 0 248 268 
6th iteraction cluster_center_new is:
63 6 112 174 204 192 104 268 
7th iteraction cluster_center_new is:
102 63 184 136 212 248 270 268 
8th iteraction cluster_center_new is:
63 198 136 184 220 248 268 270 
9th iteraction cluster_center_new is:
63 136 184 198 192 240 268 270 
10th iteraction cluster_center_new is:
63 136 184 248 198 240 268 270 
11th iteraction cluster_center_new is:
63 136 184 198 20 168 268 270 
12th iteraction cluster_center_new is:
86 63 136 220 184 198 268 270 
13th iteraction cluster_center_new is:
56 28 136 184 198 248 268 270 
14th iteraction cluster_center_new is:
30 112 136 184 260 240 268 270 
15th iteraction cluster_center_new is:
18 146 198 184 240 262 268 270 
16th iteraction cluster_center_new is:
231 156 174 198 240 262 268 270 
17th iteraction cluster_center_new is:
184 174 198 231 240 262 268 270 
18th iteraction cluster_center_new is:
28 184 198 231 40 262 268 270 
19th iteraction cluster_center_new is:
248 0 184 198 136 262 268 270 
20th iteraction cluster_center_new is:
63 174 184 198 240 262 268 270 
21th iteraction cluster_center_new is:
63 212 184 198 40 262 268 270 
22th iteraction cluster_center_new is:
231 56 184 198 248 262 268 270 
23th iteraction cluster_center_new is:
94 184 88 231 164 262 268 270 
24th iteraction cluster_center_new is:
126 154 248 174 231 262 268 270 
25th iteraction cluster_center_new is:
198 152 174 231 164 262 268 270 
26th iteraction cluster_center_new is:
156 248 174 198 231 262 268 270 
27th iteraction cluster_center_new is:
184 174 198 231 164 262 268 270 
28th iteraction cluster_center_new is:
248 174 184 198 231 262 268 270 
29th iteraction cluster_center_new is:
28 184 198 231 164 262 268 270 
30th iteraction cluster_center_new is:
30 248 174 198 231 262 268 270 
31th iteraction cluster_center_new is:
210 174 198 231 164 262 268 270 
32th iteraction cluster_center_new is:
184 174 198 92 231 262 268 270 
33th iteraction cluster_center_new is:
63 212 220 198 136 262 268 270 
34th iteraction cluster_center_new is:
63 174 198 208 192 262 268 270 
35th iteraction cluster_center_new is:
63 174 248 198 8 262 268 270 
36th iteraction cluster_center_new is:
136 63 220 144 164 172 268 270 
37th iteraction cluster_center_new is:
63 136 168 248 216 156 268 270 
38th iteraction cluster_center_new is:
63 136 208 220 30 248 268 270 
39th iteraction cluster_center_new is:
18 63 136 46 184 248 268 270 
40th iteraction cluster_center_new is:
30 92 63 220 184 248 268 270 
41th iteraction cluster_center_new is:
136 126 28 174 192 40 268 270 
42th iteraction cluster_center_new is:
30 0 156 136 174 248 268 270 
43th iteraction cluster_center_new is:
63 30 136 140 174 92 268 270 
44th iteraction cluster_center_new is:
18 126 30 136 208 212 268 270 
45th iteraction cluster_center_new is:
56 30 184 136 0 228 268 270 
46th iteraction cluster_center_new is:
63 30 112 174 184 220 268 270 
47th iteraction cluster_center_new is:
38 63 136 174 184 248 268 270 
48th iteraction cluster_center_new is:
198 63 156 174 140 240 268 270 
49th iteraction cluster_center_new is:
63 144 192 174 198 240 268 270 
50th iteraction cluster_center_new is:
63 208 212 184 198 40 100 270 
51th iteraction cluster_center_new is:
231 56 100 184 198 18 220 270 
52th iteraction cluster_center_new is:
128 56 46 174 198 248 231 270 
53th iteraction cluster_center_new is:
6 56 184 174 88 231 164 270 
54th iteraction cluster_center_new is:
253 94 122 248 174 184 136 270 
55th iteraction cluster_center_new is:
130 220 136 174 184 248 253 270 
56th iteraction cluster_center_new is:
156 136 174 184 144 248 253 270 
57th iteraction cluster_center_new is:
136 228 176 174 48 248 253 270 
58th iteraction cluster_center_new is:
128 136 174 184 0 114 112 270 
59th iteraction cluster_center_new is:
63 256 248 208 136 174 232 270 
60th iteraction cluster_center_new is:
63 198 174 46 30 248 256 270 
61th iteraction cluster_center_new is:
136 30 63 208 198 164 138 270 
62th iteraction cluster_center_new is:
18 63 136 140 248 198 46 270 
63th iteraction cluster_center_new is:
30 92 63 136 208 198 240 264 
64th iteraction cluster_center_new is:
18 142 184 192 198 208 240 264 
65th iteraction cluster_center_new is:
56 160 184 248 198 0 40 264 
66th iteraction cluster_center_new is:
63 136 102 192 174 198 32 264 
67th iteraction cluster_center_new is:
231 200 58 136 174 248 198 264 
68th iteraction cluster_center_new is:
63 192 174 198 222 8 240 264 
69th iteraction cluster_center_new is:
90 63 174 248 88 156 40 264 
70th iteraction cluster_center_new is:
231 56 134 102 220 174 164 264 
71th iteraction cluster_center_new is:
94 270 136 192 174 140 231 264 
72th iteraction cluster_center_new is:
130 136 208 174 248 231 264 270 
73th iteraction cluster_center_new is:
198 136 174 48 8 248 264 270 
74th iteraction cluster_center_new is:
253 128 192 174 198 240 264 270 
75th iteraction cluster_center_new is:
184 174 248 88 40 253 264 270 
76th iteraction cluster_center_new is:
63 28 174 140 30 136 264 270 
77th iteraction cluster_center_new is:
156 30 63 136 208 174 264 270 
78th iteraction cluster_center_new is:
18 63 136 220 174 208 264 270 
79th iteraction cluster_center_new is:
56 63 136 174 260 248 264 270 
80th iteraction cluster_center_new is:
140 63 220 174 248 262 264 270 
81th iteraction cluster_center_new is:
63 140 174 192 248 262 264 270 
82th iteraction cluster_center_new is:
63 140 174 192 248 262 264 270 
Paral KMediods++ compute the Cluster Centers Spent: 0.8762 s
```

If the `-w` parameter is 0, it will output:
```shell
Just run for test, no results output
```
before
```shell
Paral KMediods++ compute the Cluster Centers Spent: 0.8762 s
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