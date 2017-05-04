% cell_id_set={'MCF7'};  
% pert_set={'atorvastatin','vemurafenib','venlafaxine'};
% pert_type_set = {'trt_cp'};
% duration = 24;  
% concentration= 10;

% default parameters
if ~exist('sample_conditions_chd')
	sample_conditions_chd = {'cell_id', 'pert_iname', 'pert_type', 'pert_itime', 'pert_idose'};
	%sample_conditions_chd = {'CL_Name', 'SM_Name', 'SM_Pert_Type', 'SM_Time', 'SM_Dose'}; 
end

if ~exist('file_input')
	file_input='../data/modzs_n272x978.gctx';
	%file_input='../data/GSE70138_Broad_LINCS_Level2_GEX_n78980x978_2015-06-30.gct';
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

probe=1:m;            %probe index
probe=probe'; 
   
fid1 = fopen('../data/tmp', 'w');  
fid2 = fopen(file_name_cidnum, 'w');  
fid3 = fopen(file_name, 'w');

%pre-sort and extract the fit profile
count = 0;
o=ones(m,2);


if isCell && isPert && isType && isDura && isCon   %no filter conditions
	for i = 1:n
	
		count = count+1;  %count the number of fit profile
		o = [mat(:,i),probe];
		o = sortrows(o,1);		
		strprofile='';
		for j = 1:m-1   %merge the profile string
			strprofile=sprintf('%s%5d\t',strprofile,o(j,2));
		end
		strprofile=sprintf('%s%5d\n',strprofile,o(m,2));
		fprintf(fid1,'%s',strprofile);  %write out proile
		fprintf(fid2,'%10d\n',i);   %write out cid number
	end

else
	tic
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

	for i = 1:n
	
		if isDura == 0 
			if isequal(class(ds.cdesc{i,cindex(4)}),'char')
				dura_now = str2num(ds.cdesc{i,cindex(4)}(isstrprop(ds.cdesc{i,cindex(4)},'digit')));  %extract the number part
			else
				dura_now = ds.cdesc{i,cindex(4)};
			end
		end
		if isCon == 0 
			if isequal(class(ds.cdesc{i,cindex(4)}),'char')
				con_now = str2num(ds.cdesc{i,cindex(5)}(isstrprop(ds.cdesc{i,cindex(5)},'digit')));  %extract the number part
			else
				con_now = ds.cdesc{i,cindex(5)};
			end
		end
	
		if ( isCell || ismember(ds.cdesc{i,cindex(1)},cell_id_set) ) && ( isPert || ismember(ds.cdesc{i,cindex(2)}, pert_set) ) && ( isType || ismember(ds.cdesc{i,cindex(3)}, pert_type_set) ) && ( isDura || dura_now == duration ) && ( isCon || con_now == concentration )
			count = count+1;  %count the number of fit profile
			o = [mat(:,i),probe];
			o = sortrows(o,1);		
			strprofile='';
			for j = 1:m-1   %merge the profile string
				strprofile=sprintf('%s%5d\t',strprofile,o(j,2));
			end
			strprofile=sprintf('%s%5d\n',strprofile,o(m,2));
			fprintf(fid1,'%s',strprofile);  %write out proile
			fprintf(fid2,'%10d\n',i);   %write out cid number
		end
	end
end

fprintf(fid3,'%10d\t%10d\n', count,m);
toc
fclose(fid1);
fclose(fid2);
fclose(fid3);