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
ds = parse_gctx(file_input);
%ds = parse_gct(file_input);	%-------------------please change this line if the file is .gct--------------------%
mat = ds.mat;
[m,n]=size(ds.mat);

tic

probe=1:m;            %probe index
probe=probe'; 
   
fid1 = fopen('../data/tmp', 'w');  
fid2 = fopen(file_name_cidnum, 'w');  
fid3 = fopen(file_name, 'w');

%pre-sort and extract the fit profile
count = 0;
o=ones(m,2);
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
fprintf(fid3,'%10d\t%10d\n', count,m);
toc

fclose(fid1);
fclose(fid2);
fclose(fid3);