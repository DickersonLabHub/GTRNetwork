%Match Yeast Genes of C8 data
fid=fopen('SGD_features.tab');
n=1;
l=fgetl(fid);
while (l~=-1)
    s=strread(l,'%s','delimiter','\t');
    if ~isempty(s{5,1})
        Gene_bank(n,:)=[s(4,1) s(5,1) s(16,1)];
    n=n+1;
    end
        l=fgetl(fid);
end
fclose(fid);

for i=1:5814
    s=strmatch(C8_Gene(i,1),Gene_bank(:,2),'exact');
    if isempty(s)
        s=strmatch(strread(C8_Gene{i,1},'%s',1,'delimiter',' '),Gene_bank(:,2),'exact');
    end
    if ~isempty(s)
        C8_Gene(i,2)=Gene_bank(s,1);
    else
        C8_Gene(i,2)=C8_Gene(i,1);
    end
end