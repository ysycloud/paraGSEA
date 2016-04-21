<a name="ES_Matrix_ompi_p2p.doc"></a>
# ES_Matrix_ompi_p2p #

This tool Compares two gene profile sets to get an Enrichment Score matrix 
of every gene profile pairs and write out the result by MPI/OpenMP with p2p
communication.

There are several parameters should be setted in command line example.

| parameter name | Parameter function |
| -------------- | -------------------|
| -n process_num | Total number of processes |
| -ppn pernum |the number of processes in each node |
| -hostfile hostfile | list the IP or Hostname of nodes |
| thread_num | the number of threads in per process_num |
| siglen | the length of Gene Expression Signature |
| filename1 | a parsed profiles's file from pretreatment stage |
| filename2 | another parsed profiles's file from pretreatment stage |
| filename3 | output file ,distributed in every nodes ,with ES Matrix |

A sample `Shell script` file is given below that makes use of `ES_Matrix_ompi_p2p`.

```shell
mpirun -n 2 -ppn 2 -hostfile example/hostfile ./bin/ES_Matrix_ompi_p2p 4 50 "data/data_for_test.txt" "data/data_for_test.txt" "data/ES_Matrix_test"
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
loading IO and prework time in peer to peer way: 0.4629 s
Paral compute the ES_Matrix is Starting...!
Paral compute the ES_Matrix time : 2.9656 s
Writing file is Starting...!
Write Result spent: 0.0453 s
```

There is no more need of you to input anything in command line. However,
the ES_Matrix will be written to file `data/ES_Matrix_test_*.txt` in every
processes in a distributed way.