
ds = parse_gctx(file_input);
mat = ds.mat;
cid = ds.cid;
rid = ds.rid;
[m,n]=size(ds.mat);

probe=1:m;            %probe index
probe=probe'; 
   
fid1 = fopen(file_name, 'w');  %file point
fid2 = fopen(file_name_cid, 'w');  %file point
fid3 = fopen(file_name_rid, 'w');  %file point
fprintf(fid1,'%10g\t%10g\n', n,m);
tic
%pre-sort
o=ones(m,2);
for i = 1:n
    o = [mat(:,i),probe];
    o = sortrows(o,1);
    for j = 1:m-1
         fprintf(fid1,'%5g\t',o(j,2));
    end
    fprintf(fid1,'%5g\n',o(m,2));
	fprintf(fid2,'%s\n',char(cid(i)));
end
for i = 1:m
	fprintf(fid3,'%d->%s\n',i,char(rid(i)));
end
toc
fclose(fid1);
fclose(fid2);
fclose(fid3);

