#/bin/bash

# execute matlab script to parallel parse the data  
#matlab -nodesktop -nosplash -nojvm -r " file_input='../data/modzs_n272x978.gctx';  file_name='../data/data_for_test.txt';  file_name_cid='../data/data_for_test_cidnum.txt';  cores = 2;  cell_id_set={'A549','MCF7','A375','A673','AGS'};  pert_set={'BRD-A51714012-001-04-9','BRD-A51714012-001-03-1','BRD-U88459701-000-01-8'};  duration = '6H';  concentration='10';  paraPreGSEA;  quit;"

# execute matlab script to parse the data by acquiescent parameters
matlab -nodesktop -nosplash -nojvm -r "paraPreGSEA; quit;"

cat ../data/data_for_test.txt_* >> ../data/data_for_test.txt
cat ../data/data_for_test_cidnum.txt_* >> ../data/data_for_test_cidnum.txt
rm -f ../data/data_for_test.txt_* ../data/data_for_test_cidnum.txt_*

