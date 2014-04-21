function [ BindingSite ] = load_bindingsites( filename )
% Read regulonDB binding site information data file
% by Fu, Yao(Al)

fid=fopen(filename);

N=0;
l=fgetl(fid);
if size(l,1)==0
    inhead=1;
elseif l(1,1)=='#' 
     inhead=1;
else inhead=0;
end
while (inhead) %skip header
    if size(l,1)==0
        l=fgetl(fid);
    elseif ~isempty(strfind(l,'#')) || ~isempty(strfind(l,'Release'))
        l=fgetl(fid);
    else
        inhead=0;
    end
end

if l~=-1
    infile=1;
    s=strread(l,'%s','delimiter','\t');
    N=N+1;
    fprintf('%s%i\n','Loading operon number ',N);
    BindingSite(N,1)=s(2,1);
    BindingSite(N,2)=s(8,1);
else
    infile=0;
    fprintf('%s\n','File error!');
end

while (infile)
    l=fgetl(fid);
    if l~=-1
       s=strread(l,'%s','delimiter','\t');
       N=N+1;
       fprintf('%s%i\n','Loading operon number ',N);
       BindingSite(N,1)=s(2,1);
       BindingSite(N,2)=s(8,1);
    else
        infile=0;
    end
end
fclose(fid);
end


