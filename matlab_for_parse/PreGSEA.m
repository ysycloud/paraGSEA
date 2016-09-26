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

isDura = 0;
if ~exist('duration') %not set, duration will not be consider
	isDura = 1;
end

isCon = 0;
if ~exist('concentration') %not set, concentration will not be consider
	isCon = 1;
end

ds = parse_gctx(file_input);
mat = ds.mat;
cid = ds.cid;
[m,n]=size(ds.mat);

probe=1:m;            %probe index
probe=probe'; 
   
fid1 = fopen('../data/tmp', 'w');  
fid2 = fopen(file_name_cidnum, 'w');  
fid3 = fopen(file_name, 'w');

tic
count = 0;
cid_tmp = regexp(cid,':','split');
%pre-sort and extract the fit profile
o=ones(m,2);
for i = 1:n
	cid_info{i} = regexp(cid_tmp{i},'_','split');
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