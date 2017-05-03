if ~exist('file_input')
	%file_input='../data/modzs_n272x978.gctx';
	%file_input='../data/GSE70138_Broad_LINCS_Level2_GEX_n78980x978_2015-06-30.gct';
	%file_input='../../Lincs/GSE92742_Broad_LINCS_Level2_GEX_delta_n49216x978.gctx';
	file_input='../../Lincs/GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1278882x978.gctx';
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

ds = parse_gctx(file_input);
%ds = parse_gct(file_input);	%-------------------please change this line if the file is .gct--------------------%
mat = ds.mat;
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
dura_now = zeros(1,cores);
con_now = zeros(1,cores);

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
	
	%pre-sort and extract the fit profile in local thread
	for j = 1:local_n
		lcount(i)=lcount(i)+1;	 %count the number of fit profile in local thread
		o = [mat(:,begin+j),probe];
		o = sortrows(o,1);
		for k = 1:m-1	%write out the profile
			fprintf(fid1(i),'%5d\t',o(k,2));
		end
		fprintf(fid1(i),'%5d\n',o(m,2));
		fprintf(fid2(i),'%10d\n',begin+j);	 %write out cid number
	end
	fclose(fid1(i));
	fclose(fid2(i));
end

count = sum(lcount);
fid = fopen(file_name, 'w');  %file point
fprintf(fid,'%10d\t%10d\n', count,m);
fclose(fid);
toc

matlabpool close;