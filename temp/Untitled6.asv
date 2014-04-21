fid = fopen('ibtest.txt','w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','#NAMES','IB1','IB2','IB3','IB4','IB5','IB6','IB7','IB8','IB9','IB10');
for i=1:size(R,1)
    fprintf(fid,'%s',prob{i,2});
    for j=5:14
        if R(i,j)==0
            fprintf(fid,'\t');
        else
        fprintf(fid,'\t%f',R(i,j));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('controltest.txt','w');
for i=1:size(R,1)
    fprintf(fid,'%s',prob{i,2});
    for j=1:4
        if R(i,j)==0
            fprintf(fid,'\t');
        else
        fprintf(fid,'\t%f',R(i,j));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

fid=fopen('IBpreprocessing.out');
fidn=fopen('IB.txt','w');
l=fgetl(fid);
while (l~=-1)
    if isempty(strmatch('#',l))
       fprintf(fidn,'%s\n',l);
    end
    l=fgetl(fid);
end
fclose(fid);
fclose(fidn);

n=1;
for i=1:size(CID,1)
    s=strmatch(CID(i,1),IBid(:,1),'exact');
    if ~isempty(s)
           Expression_Butanol.Gene_id(n,1)=CID(i,1);
           Expression_Butanol.R(n,1:4)=RC(i,1:4);
           Expression_Butanol.R(n,5:14)=RIB(s,1:10);
           n=n+1;
       
    end
end