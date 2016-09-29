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

isType = 0;
if ~exist('pert_type_set') %not set, perturbation type will not be consider
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
%ds = parse_gct(file_input);	%-------------------please change this line if the file is .gct--------------------%
mat = ds.mat;
[m,n]=size(ds.mat);

%check samples field number
[cisin,cindex] = ismember(sample_conditions_chd, ds.chd);

[~,fn] = size(cisin);
if fn~=5   %only judge 5 fields
	disp('[ Field Error ]sample conditions field number error, please make sure before generating a reference!');
	clear;
	return;
end

if cisin(1)==0
	disp('cell line field name error, we will ignore cell_id_set');
	isCell = 1;
end
if cisin(2)==0
	disp('perturbation field name error, we will ignore pert_set');
	isPert = 1;
end
if cisin(3)==0
	disp('perturbation type field name error, we will ignore pert_type_set');
	isType = 1;
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
	
		if isDura == 0 
			dura_now(i) = str2num(ds.cdesc{begin+j,cindex(4)}(isstrprop(ds.cdesc{begin+j,cindex(4)},'digit')));  %extract the number part
		end
		if isCon == 0 
			con_now(i) = str2num(ds.cdesc{begin+j,cindex(5)}(isstrprop(ds.cdesc{begin+j,cindex(5)},'digit')));  %extract the number part
		end
		
		if ( isCell || ismember(ds.cdesc{begin+j,cindex(1)},cell_id_set) ) && ( isPert || ismember(ds.cdesc{begin+j,cindex(2)}, pert_set) ) && ( isType || ismember(ds.cdesc{begin+j,cindex(3)}, pert_type_set) ) && ( isDura || dura_now(i) == duration ) && ( isCon || con_now(i) == concentration )
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