function [ Operon ] = load_operon( filename )
% Read Operon information from regulonDB data file
% by Fu, Yao(Al)
% Modified on 12/15/2010 to be able to read RegulonDB 7.0

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
    Operon.id(N,1)=s(1,1);
    Operon.n(N,1)=str2num(s{2,1});
    g=strread(s{4,1},'%s','delimiter',',');
    for i=1:Operon.n(N,1)
        Operon.gene(N,i)=strread(g{i,1},'%s',1,'whitespace','|');
    end
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
       Operon.id(N,1)=s(1,1);
        Operon.n(N,1)=str2num(s{2,1});
        g=strread(s{4,1},'%s','delimiter',',');
        for i=1:Operon.n(N,1)
            Operon.gene(N,i)=strread(g{i,1},'%s',1,'whitespace','|');
        end            
    else
        infile=0;
    end
end
fclose(fid);
end

