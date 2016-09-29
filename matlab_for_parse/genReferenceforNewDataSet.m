if ~exist('gene_symbol_rhd')
	gene_symbol_rhd = 'pr_gene_symbol';
end

if ~exist('sample_conditions_chd')
	sample_conditions_chd = {'cell_id', 'pert_iname', 'pert_type', 'pert_itime', 'pert_idose'}; 
end

if ~exist('datasource')
	datasource='../data/modzs_n272x978.gctx';
end

ds= parse_gctx(datasource);
%ds= parse_gct(datasource);
[m,n]=size(ds.mat);

%write out GeneSymbols corresponding to the rid in new dataset 
[risin,rindex] = ismember(gene_symbol_rhd, ds.rhd);
fid1 = fopen('../data/Reference/Gene_List.txt', 'w');
for i=1:m
	fprintf(fid1,'%s\n',ds.rdesc{i,rindex});
end
fclose(fid1);

%write out cid/condition_info and line offset file in new dataset
[cisin,cindex] = ismember(sample_conditions_chd, ds.chd);
fid2 = fopen('../data/Reference/Samples_Condition.txt', 'w');
fid3 = fopen('../data/Reference/Samples_RowByteOffset.txt', 'w');
offset = 0;
for i=1:n
	fprintf(fid3,'%10d\t',offset);
	if cisin(1)==0
		cell_line = 'field error';
	else
		cell_line = ds.cdesc{i,cindex(1)};
	end
	if cisin(2)==0
		perturbation = 'field error';
	else
		perturbation = ds.cdesc{i,cindex(2)};
	end
	if cisin(3)==0
		perturbation_type = 'field error';
	else
		perturbation_type = ds.cdesc{i,cindex(3)};
	end
	if cisin(4)==0
		duration = 'field error';
	else
		duration = ds.cdesc{i,cindex(4)};
	end
	if cisin(5)==0
		concentration = 'field error';
	else
		concentration = ds.cdesc{i,cindex(5)};
	end
	fprintf(fid2,'cell_line: %10s;    perturbation: %15s;    perturbation type: %10s;    duration: %10s;    concentration: %10s\n', cell_line, perturbation, perturbation_type, duration, concentration );
	offset = ftell(fid2);
end
fclose(fid2);
fclose(fid3);
