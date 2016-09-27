#/bin/bash

# execute matlab script to generate reference directory
#matlab -nodesktop -nosplash -nojvm -r " datasource='../data/modzs_n272x978.gctx'; genReferenceforNewDataSet; quit; "

# execute matlab script to generate reference directory by acquiescent parameters
matlab -nodesktop -nosplash -nojvm -r " genReferenceforNewDataSet; quit;"