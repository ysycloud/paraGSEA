
ds = parse_gctx(file_input);
mat = ds.mat;
cid = ds.cid;
rid = ds.rid;
[m,n]=size(ds.mat);

probe=1:m;            %probe index
probe=probe'; 

%open的matlabpool
if matlabpool('size')<=0 
    matlabpool('open','local',cores);
else
    disp('Already initialized');
end

fid = fopen(file_name, 'w');  %file point
fid3 = fopen(file_name_rid, 'w');  %file point
fprintf(fid,'%10g\t%10g\n', n,m);
tic
%pre-sort
local = floor(n/cores);
leave = mod(n,cores);
fid1 = ones(1,cores);
fid2 = ones(1,cores); 

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
	fid2(i) = fopen(sprintf('%s_%d',file_name_cid,i), 'w');
	for j = 1:local_n
		o = [mat(:,begin+j),probe];
		o = sortrows(o,1);
		for k = 1:m-1
			fprintf(fid1(i),'%5g\t',o(k,2));
		end
		fprintf(fid1(i),'%5g\n',o(m,2));
		fprintf(fid2(i),'%s\n',char(cid(begin+j)));
	end
	fclose(fid1(i));
	fclose(fid2(i));
end

for i = 1:m
	fprintf(fid3,'%d->%s\n',i,char(rid(i)));
end
toc

fclose(fid);
fclose(fid3);
matlabpool close;