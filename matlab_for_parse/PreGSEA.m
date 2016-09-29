% default parameters
if ~exist('sample_conditions_chd')
	sample_conditions_chd = {'cell_id', 'pert_iname', 'pert_type', 'pert_itime', 'pert_idose'}; 
end

if ~exist('file_input')
	file_input='../data/modzs_n272x978.gctx';
end

if ~exist('file_name')
	file_name='../data/data_for_test.txt';
end

if ~exist('file_name_cidnum')
	file_name_cidnum='../data/data_for_test_cidnum.txt';
end

isCell = 0;
if ~exist('cell_id_set') %not set, cell will not be consider
	isCell = 1;
end

isPert = 0;
if ~exist('pert_set') %not set, perturbation will not be consider
	isPert = 1;
end

isType = 0;
if ~exist('pert_type') %not set, perturbation will not be consider
	isType = 1;
end

isDura = 0;
if ~exist('duration') %not set, duration will not be consider
	isDura = 1;
end

isCon = 0;
if ~exist('concentration') %not set, concentration will not be consider
	isCon = 1;
end

ds = parse_gctx(file_input);
%ds = parse_gct(file_input);
mat = ds.mat;
[m,n]=size(ds.mat);

tic
%check the field name
[cisin,cindex] = ismember(sample_conditions_chd, ds.chd);
if cisin(1)==0
	disp('cell line field name error, we will ignore cell_id_set');
	isCell = 1;
end
if cisin(2)==0
	disp('perturbation field name error, we will ignore pert_set');
	isCell = 1;
end
if cisin(3)==0
	disp('perturbation type field name error, we will ignore pert_type');
	isPert = 1;
end
if cisin(4)==0
	disp('duration field name error, we will ignore duration');
	isDura = 1;
end
if cisin(5)==0
	disp('concentration field name error, we will ignore concentration');
	isCon = 1;
end

probe=1:m;            %probe index
probe=probe'; 
   
fid1 = fopen('../data/tmp', 'w');  
fid2 = fopen(file_name_cidnum, 'w');  
fid3 = fopen(file_name, 'w');

%pre-sort and extract the fit profile
count = 0;
o=ones(m,2);
for i = 1:n
	
	
	if ( isCell || ismember(cid_info{i}{1}{2},cell_id_set) ) && ( isPert || ismember(cid_info{i}{2}{1},pert_set) ) && ( isDura || isequal(cid_info{i}{1}{3},duration) ) && ( isCon || isequal(cid_info{i}{3}{1},concentration) )
		count = count+1;  %count the number of fit profile
		o = [mat(:,i),probe];
		o = sortrows(o,1);
		for j = 1:m-1   %write out the profile
			fprintf(fid1,'%5g\t',o(j,2));
		end
		fprintf(fid1,'%5g\n',o(m,2));
		fprintf(fid2,'%10d\n',i);   %write out cid number
	end
end
fprintf(fid3,'%10g\t%10g\n', count,m);
toc

fclose(fid1);
fclose(fid2);
fclose(fid3);