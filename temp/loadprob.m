function [ prob_data ] = loadprob( filename )
fid=fopen(filename);
l=fgetl(fid);
n=1;
while(l~=-1);
    [a b c]=strread(l,'%s%*s%*s%*s%*s%s%s',7,'delimiter','\t');
    if size(a,1)==1
        prob_data(n,[1 2 3])=[a b c];
        n=n+1;
    end
    l=fgetl(fid);      
end
fclose(fid);
end