function [ TF_gene] = load_TF_gene_Yeast( filename,sigma_file )
%  Read RegulonDB TF-gene interaction File
%   C: connetctivity matrix
%   AI: list All information
%       column 1: TF name, 2:Regulated Gene,3:TF gene, 4:Regulation Type 5:Detail
%   CD: Connectivities with Direction
%   By Fu, Yao (Al)

fid=fopen(filename);

G=0;
T=0;
N=0;
l=fgetl(fid);
if l(1,1)=='#'
    inhead=1;
else inhead=0;
end
while (inhead) %skip header
    if l(1,1)=='#'
        l=fgetl(fid);
    else
        inhead=0;
    end
end

if l~=-1
    infile=1;
    s=strread(l,'%s','delimiter','\t');
    G=G+1;
    T=T+1;
    N=N+1;
    fprintf('%s%i\n','Loading TF gene link number ',N);
    TF_gene.Gene_id(G,1)=s(2,1);
    TF_gene.TF_id(T,1)=s(1,1)';
    TF_gene.AI(N,1:2)=s([1 2],1)';
    TF_gene.C(N,G)=1;
%    TF_gene.CD(N,G)=s(6,1);
else
    infile=0;
    fprintf('%s\n','File error!');
end

while (infile)
    l=fgetl(fid);
    if l~=-1
       s=strread(l,'%s','delimiter','\t');
        N=N+1;
        TF_gene.AI(N,1:2)=s([1 2],1)';
        fprintf('%s%i\n','Loading TF gene link number ',N);
       g=strmatch(s(2,1),TF_gene.Gene_id(:,1),'exact');
       if isempty(g)
           G=G+1;
           g=G;
           TF_gene.Gene_id(G,1)=s(2,1);
       end
       t=strmatch(s(1,1),TF_gene.TF_id(:,1),'exact');
       if isempty(t)
           T=T+1;
           t=T;
           TF_gene.TF_id(T,1)=s(1,1)';
       end
       TF_gene.C(t,g)=1;
%       TF_gene.CD(t,g)=s(6,1);             
    else
        infile=0;
    end
end
fclose(fid);
TF_gene.ref=TF_gene.C;
%TF_gene.refC=TF_gene.CD;
M=N;

if exist('sigma_file') %Loading Sigma factor gene links
    fid = fopen(sigma_file);
    l=fgetl(fid);
    if l(1,1)=='#'
    inhead=1;
    else inhead=0;
    end
    while (inhead) %skip header
        if l(1,1)=='#'
          l=fgetl(fid);
        else
            inhead=0;
        end
    end
    
    if l~=-1
        infile=1;
       s=strread(l,'%s','delimiter','\t');
        N=N+1;
        TF_gene.AI(N,[1 2 4])=s([1 2 3],1)';
        fprintf('%s%i\n','Loading Sigma gene link number ',N-M);
       g=strmatch(s(2,1),TF_gene.Gene_id(:,1),'exact');
       if isempty(g)
           G=G+1;
           g=G;
           TF_gene.Gene_id(G,1:2)=[s(2,1) s(4,1)];
       end
       t=strmatch(s(1,1),TF_gene.TF_id(:,1),'exact');
       if isempty(t)
           T=T+1;
           t=T;
           TF_gene.TF_id(T,1)=s(1,1)';
       end
       TF_gene.C(t,g)=1;
       TF_gene.CD(t,g)=s(3,1);             
    else
        infile=0;
        fprintf('%s\n','File error!');
    end
    
    while (infile)
    l=fgetl(fid);
    if l~=-1
       s=strread(l,'%s','delimiter','\t');
        N=N+1;
        TF_gene.AI(N,[1 2 4])=s([1 2 3],1)';
        fprintf('%s%i\n','Loading Sigma gene link number ',N-M);
       g=strmatch(s(2,1),TF_gene.Gene_id(:,1),'exact');
       if isempty(g)
           G=G+1;
           g=G;
           TF_gene.Gene_id(G,1:2)=[s(2,1) s(4,1)];
       end
       t=strmatch(s(1,1),TF_gene.TF_id(:,1),'exact');
       if isempty(t)
           T=T+1;
           t=T;
           TF_gene.TF_id(T,1)=s(1,1)';
       end
       TF_gene.C(t,g)=1;
       TF_gene.CD(t,g)=s(3,1);             
    else
        infile=0;
    end
    end
    fclose(fid);
    TF_gene.ref=TF_gene.C;
    TF_gene.refC=TF_gene.CD;

end

end


