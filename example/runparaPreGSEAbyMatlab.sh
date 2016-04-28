#/bin/bash

# execute matlab script to parallel parse the data  
matlab -nodesktop -nosplash -nojvm -r "file_input='../data/modzs_n272x978.gctx'; file_name='../data/data_for_test.txt'; file_name_cid='../data/data_for_test_cid.txt'; file_name_rid='../data/data_for_test_rid.txt'; cores = 2; paraPreGSEA; quit;"

cat ../data/data_for_test.txt_* >> ../data/data_for_test.txt
cat ../data/data_for_test_cid.txt_* >> ../data/data_for_test_cid.txt
rm -f ../data/data_for_test.txt_* ../data/data_for_test_cid.txt_*

