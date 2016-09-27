if ~exist('datasource')
	datasource='../data/modzs_n272x978.gctx';
end

load annot.mat;
annot_info = regexp(annot,'\t','split');
clear annot
%extract GeneSymbol
[at_m,~]=size(annot_info);    
for i=1:at_m
	annot_Symbol{i} = annot_info{i}{3};
	annot_rid{i} = annot_info{i}{1};
end

ds= parse_gctx(datasource);
[m,n]=size(ds.mat);
ds_rid = ds.rid;
ds_cid = ds.cid;

%write out GeneSymbols corresponding to the rid in new dataset 
[local,index] = ismember(annot_rid,ds_rid);
o=ones(at_m,2);
probe=1:at_m; 
o=[index',probe'];
o= sortrows(o,1);
fid1 = fopen('../data/Reference/Gene_List.txt', 'w');
for i=1:m
	j = at_m-m+i;
	geneSymbol(i)=annot_Symbol(o(j,2));
	fprintf(fid1,'%s\n',geneSymbol{i});
end
fclose(fid1);

%write out cid/condition_info and line offset file in new dataset 
fid2 = fopen('../data/Reference/Samples_Condition.txt', 'w');
fid3 = fopen('../data/Reference/Samples_RowByteOffset.txt', 'w');
offset = 0;
cid_tmp = regexp(ds_cid,':','split');
for i=1:n
	cid_info{i} = regexp(cid_tmp{i},'_','split');
	fprintf(fid3,'%10d\t',offset);
	fprintf(fid2,'cid: %s;\tcell_line: %s;\tperturbation: %s;\tduration: %s;\tconcentration: %sum\n',ds_cid{i},cid_info{i}{1}{2},cid_info{i}{2}{1},cid_info{i}{1}{3},cid_info{i}{3}{1});
	offset = ftell(fid2);
end
fclose(fid2);
fclose(fid3);