if ~exist('file_input')
	file_input='../data/modzs_n272x978.gctx';
end

if ~exist('file_name')
	file_name='../data/data_for_test.txt';
end

if ~exist('file_name_cidnum')
	file_name_cidnum='../data/data_for_test_cidnum.txt';
end

if ~exist('cores')
	cores=2;
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

%open的matlabpool
if matlabpool('size')<=0 
    matlabpool('open','local',cores);
else
    disp('Already initialized');
end


tic
%pre-sort
local = floor(n/cores);
leave = mod(n,cores);
fid1 = ones(1,cores);
fid2 = ones(1,cores);
lcount = zeros(1,cores);
 
cid_tmp = regexp(cid,':','split');
for i=1:n
	cid_info{i} = regexp(cid_tmp{i},'_','split');
end

parfor i=1:cores
	o=ones(m,2);
	local_n = local;
	begin = 1;
	if i <= leave
		local_n = local_n+1;
		begin = (i-1)*local_n;
	else
		begin = (i-1)*local_n+leave;
	end
	fid1(i) = fopen(sprintf('%s_%d',file_name,i), 'w');
	fid2(i) = fopen(sprintf('%s_%d',file_name_cidnum,i), 'w');
	count = 0;
	
	%pre-sort and extract the fit profile in local thread
	for j = 1:local_n
		if ( isCell || ismember(cid_info{begin+j}{1}{2},cell_id_set) ) && ( isPert || ismember(cid_info{begin+j}{2}{1},pert_set) ) && ( isDura || isequal(cid_info{begin+j}{1}{3},duration) ) && ( isCon || isequal(cid_info{begin+j}{3}{1},concentration) )
			lcount(i)=lcount(i)+1;	 %count the number of fit profile in local thread
			o = [mat(:,begin+j),probe];
			o = sortrows(o,1);
			for k = 1:m-1	%write out the profile
				fprintf(fid1(i),'%5g\t',o(k,2));
			end
			fprintf(fid1(i),'%5g\n',o(m,2));
			fprintf(fid2(i),'%10d\n',begin+j);	 %write out cid number
		end
	end
	fclose(fid1(i));
	fclose(fid2(i));
end

count = sum(lcount);
fid = fopen(file_name, 'w');  %file point
fprintf(fid,'%10g\t%10g\n', count,m);
fclose(fid);
toc

matlabpool close;