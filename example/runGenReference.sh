#/bin/bash

# execute matlab script to generate reference directory
#matlab -nodesktop -nosplash -nojvm -r " datasource='../data/modzs_n272x978.gctx'; gene_symbol_rhd = 'pr_gene_symbol'; sample_conditions_chd = {'cell_id', 'pert_iname', 'pert_type', 'pert_itime', 'pert_idose'}; genReferenceforNewDataSet; quit; "

# execute matlab script to generate reference directory by acquiescent parameters
matlab -nodesktop -nosplash -nojvm -r " genReferenceforNewDataSet; quit;"