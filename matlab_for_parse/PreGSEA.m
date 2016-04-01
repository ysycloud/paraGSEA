
ds = parse_gctx(file_input);

mat = ds.mat;
cid = ds.cid;
[m,n]=size(ds.mat);

probe=1:m;            %probe index
probe=probe'; 
   
fid1 = fopen(file_name, 'w');  %file point1
fid2 = fopen(file_name_cid, 'w');  %file point1
fprintf(fid1,'%10g\t%10g\n', n,m);
tic
%pre-sort
o=ones(m,2);
for i = 1:n
    o = [mat(:,i),probe];
    o = sortrows(o,1);
    mat(:,i)=o(:,2);
    for j = 1:m-1
         fprintf(fid1,'%5g\t',o(j,2));
    end
    fprintf(fid1,'%5g\n',o(m,2));
	fprintf(fid2,'%s\n',char(cid(i)));
end
toc
fclose(fid1);
fclose(fid2);

