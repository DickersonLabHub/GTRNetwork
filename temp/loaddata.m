function [ exp_data ] = loaddata( filename )
fid=fopen(filename);
l=fgetl(fid);
n=1;
while(l~=-1);
    s=strread(l,'%s',2);
    if ~strcmp('null',s{2,1})
        exp_data(n,[1 2])=strread(l,'%f',2);
        n=n+1;
    end
    l=fgetl(fid);      
end
fclose(fid);
end