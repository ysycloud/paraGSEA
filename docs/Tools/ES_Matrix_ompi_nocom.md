<a name="ES_Matrix_ompi_nocom.doc"></a>
# ES_Matrix_ompi_nocom #

This tool Compares two gene profile sets to get an Enrichment Score matrix 
of every gene profile pairs and write out the result by MPI/OpenMP with no
communication when distribute the input file.

There are several parameters should be setted in command line example.

| parameter name | Parameter function |
| -------------- | -------------------|
| -n process_num | Total number of processes |
| -ppn pernum |the number of processes in each node |
| -hostfile hostfile | list the IP or Hostname of nodes |
| -t --thread | define the maximum number of parallel threads |
| -l--siglen | define the length of Gene Expression Signature |
| -a --loadtime | define the load time of dataset2 |
| -p --proportion | define the proportion of dataset be used |
| -w --write | decide whether output the results  |
| -1 --input1 | a parsed profiles's file from pretreatment stage |
| -2 --input2 | another parsed profiles's file from pretreatment stage |
| -o --output | output file ,distributed in every nodes ,with ES Matrix |

A sample `Shell script` file is given below that makes use of `ES_Matrix_ompi_nocom`.

```shell
mpirun -n 2 -ppn 2 -hostfile example/hostfile ES_Matrix_ompi_nocom -t 4 -l 50 -a 1 -p 1 -w 1 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix_test
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
Available Memory:      522306848 KB
Needed Memory:      942 KB
All Needed Memory:      1884 KB
phase 1  --> loading IO and prework time in no communication way: 0.2017 s
phase 1  --> Paral compute the ES_Matrix is Starting...!
phase 1  --> Paral compute the ES_Matrix time : 2.7954 s
Writing file is Starting...!
Write Result spent: 0.0390 s
```

If Available Memory is less than Needed Memory, it may produce the following output:
```shell
available memory is not enough to store all results, recommend to use more than 10 nodes!!!
```
the number of node recommended to use is gotten by valid calculation according to the machine environment and the scale of dataset.
Using it can always solve the memory shortage problem.

If the `-w` parameter is 0, even the Available Memory is not enough, it will still carry out the calculation process and may produce the following output:
```shell
available memory is not enough to store all results, recommend to use more than 10 nodes!!!
because we are just testing without writing, we will continue!!!
phase 1  --> loading IO and prework time in collective comunication: 0.4318 s
phase 1  --> Paral compute the ES_Matrix is Starting...!
phase 1  --> Paral compute the ES_Matrix time : 2.8079 s
Just run for test, no results output.
```

There is no more need of you to input anything in command line. However,
the ES_Matrix will be written to file `data/ES_Matrix_test_*.txt` in every
processes in a distributed way if the `-w` parameter is 1.