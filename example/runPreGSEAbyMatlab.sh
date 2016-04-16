#!/bin/bash

# execute matlab script to parse the data  
matlab -nodesktop -nosplash -nojvm -r "file_input='../data/modzs_n272x978.gctx'; file_name='../data/data_for_test.txt'; file_name_cid='../data/data_for_test_cid.txt'; PreGSEA; quit;"


