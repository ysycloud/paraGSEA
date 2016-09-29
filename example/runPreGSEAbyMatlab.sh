#/bin/bash

# execute matlab script to parse the data  
#matlab -nodesktop -nosplash -nojvm -r " file_input='../data/modzs_n272x978.gctx';  file_name='../data/data_for_test.txt';  file_name_cidnum='../data/data_for_test_cidnum.txt';  cell_id_set={'A549','MCF7','A375','A673','AGS'};  pert_set={'atorvastatin','vemurafenib','venlafaxine'}; trt_cp = 'trt_cp'; duration = 6 ;  concentration= 10;  PreGSEA;  quit;"

# execute matlab script to parse the data by acquiescent parameters
matlab -nodesktop -nosplash -nojvm -r "PreGSEA; quit;"

cat ../data/tmp >> ../data/data_for_test.txt
rm -f ../data/tmp