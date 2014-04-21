%Load Pathway information
fid = fopen('geneCommonNames.txt');
infile=1;
n=1;
while (infile)
    l=fgetl(fid);
    if l==-1
        infile=0;
    else
        EC2C(n,:)=strread(l,'%s',2)';
        n=n+1;
    end
end
fclose(fid);

fid = fopen('pathwayCommonNames.txt');
infile=1;
n=1;
while (infile)
    l=fgetl(fid);
    if l==-1
        infile=0;
    else
        PW2C(n,:)=strread(l,'%s',2,'delimiter','\t')';
        n=n+1;
    end
end
fclose(fid);

fid = fopen('pathwaysOfGenes.txt');
infile=1;
n=1;
while (infile)
    l=fgetl(fid);
    if l==-1
        infile=0;
    elseif size(strread(l,'%s','delimiter','\t'),1)>1
        PW2G(n,1:size(strread(l,'%s','delimiter','\t'),1))=strread(l,'%s',size(strread(l,'%s','delimiter','\t'),1),'delimiter','\t')';
        n=n+1;
    end
end
fclose(fid);

for i=1:970
    g=strmatch(PW2G(i,1),EC2C(:,1),'exact');
    PW2Gcommon(i,1)=EC2C(g,2);
    for j=2:16
        if ~isempty(PW2G{i,j})
            PW2Gcommon(i,j)=PW2C(strmatch(PW2G(i,j),PW2C(:,1),'exact'),2);
        end
    end
end

%Find Key TFs

C8=KeyTFs(ExpressionC8.E(:,1:2),TFA600.TFA(:,1:2),ExpressionC8.E(:,3:5),TFA600.TFA(:,3:5),AllLinks,ExpressionC8.Gene_id,TFA600.tf_id,0.1,0.1);
end

%Map target genes to TFs
AllLinks=Regulatory.AI(:,[1 2 4]);
n=1;
for i=1:size(AllLinks,1)
if ~isempty(strmatch(AllLinks(i,1),C8.Activated_TF(:,1),'exact'))
LowPHList(n,:)=AllLinks(i,:);
n=n+1;
elseif ~isempty(strmatch(AllLinks(i,2),C8.Activated_TF(:,1),'exact'))
LowPHList(n,:)=AllLinks(i,:);
n=n+1;
elseif ~isempty(strmatch(AllLinks(i,1),C8.Silented_TF(:,1),'exact'))
LowPHList(n,:)=AllLinks(i,:);
n=n+1;
elseif ~isempty(strmatch(AllLinks(i,2),C8.Silented_TF(:,1),'exact'))
LowPHList(n,:)=AllLinks(i,:);
n=n+1;
end
end
C8List=LowPHList;

for i=1:size(C8List,1)
if ~isempty(strmatch(C8List(i,2),C8.UpGenes,'exact'))
C8List{i,5}='+';
end
if ~isempty(strmatch(C8List(i,2),C8.DownGenes,'exact'))
C8List{i,5}='-';
end
if ~isempty(strmatch(C8List(i,1),C8.Activated_TF,'exact'))
C8List{i,6}='+';
end
if ~isempty(strmatch(C8List(i,1),C8.Silented_TF,'exact'))
C8List{i,6}='-';
end
end

% Map pathway to genes
for i=1:size(C8List,1)
p=strmatch(C8List(i,2),PW2Gcommon(:,1),'exact');
if ~isempty(p)
C8List(i,7:21)=PW2Gcommon(p,2:16);
end
end

% Export to txt file
fid = fopen('C8.txt','w');
for i=1:1242
    for j=1:21
fprintf(fid,'%s\t',C8List{i,j});
    end
    fprintf(fid,'\n');
end
fclose(fid)


% Summarize Pathway information
for i=1:size(C8List,1)
    if isempty(C8List(i,3))
        PredictDC8(i,1)=0;
    elseif strcmp(C8List{i,3},'+') && strcmp(C8List{i,6},'+')
        PredictDC8(i,1)=1;
    elseif strcmp(C8List{i,3},'+') && strcmp(C8List{i,6},'-')
        PredictDC8(i,1)=-1;
    elseif strcmp(C8List{i,3},'-') && strcmp(C8List{i,6},'+')
        PredictDC8(i,1)=-1;
    elseif strcmp(C8List{i,3},'-') && strcmp(C8List{i,6},'-')
        PredictDC8(i,1)=1;
    else
        PredictDC8(i,1)=0;
    end
    if isempty(C8List(i,5))
        PredictDC8(i,2)=0;
    elseif C8List{i,5}=='+'
        PredictDC8(i,2)=1;
    elseif C8List{i,5}=='-'
        PredictDC8(i,2)=-1;
    end
end
        

n=1;
for i=1:size(C8List,1)
    for j=7:21
        if ~isempty(C8List{i,j})
            if n==1
                Pathway(n,1)=C8List(i,j);
                Pathway{n,2}=1;
                if PredictDC8(i,1)==1
                    Pathway{n,3}=1;
                    Pathway{n,4}=0;
                elseif PredictDC8(i,1)==-1
                    Pathway{n,3}=0;
                    Pathway{n,4}=1;
                else
                    Pathway{n,3}=0;
                    Pathway{n,4}=0;
                end
                if PredictDC8(i,2)==1
                    Pathway{n,5}=1;
                    Pathway{n,6}=0;
                elseif PredictDC8(i,2)==-1
                    Pathway{n,5}=0;
                    Pathway{n,6}=1;
                else
                    Pathway{n,5}=0;
                    Pathway{n,6}=0;
                end
                Pathway(n,7)=C8List(i,1);
                n=n+1;
            else
                s=strmatch(C8List(i,j),Pathway(:,1),'exact');
                if isempty(s)
                Pathway(n,1)=C8List(i,j);
                Pathway{n,2}=1;
                if PredictDC8(i,1)==1
                    Pathway{n,3}=1;
                    Pathway{n,4}=0;
                elseif PredictDC8(i,1)==-1
                    Pathway{n,3}=0;
                    Pathway{n,4}=1;
                else
                    Pathway{n,3}=0;
                    Pathway{n,4}=0;
                end
                if PredictDC8(i,2)==1
                    Pathway{n,5}=1;
                    Pathway{n,6}=0;
                elseif PredictDC8(i,2)==-1
                    Pathway{n,5}=0;
                    Pathway{n,6}=1;
                else
                    Pathway{n,5}=0;
                    Pathway{n,6}=0;
                end
                Pathway(n,7)=C8List(i,1);
                n=n+1;
                else
                    Pathway{s,2}=Pathway{s,2}+1;
                if PredictDC8(i,1)==1
                    Pathway{s,3}=Pathway{s,3}+1;
                elseif PredictDC8(i,1)==-1
                    Pathway{s,4}=Pathway{s,4}+1;
                end
                if PredictDC8(i,2)==1
                    Pathway{s,5}=Pathway{s,5}+1;
                elseif PredictDC8(i,2)==-1
                    Pathway{s,6}=Pathway{s,6}+1;
                end
                    if isempty(strfind(Pathway{s,7},C8List{i,1}))
                        Pathway{s,7}=[Pathway{s,7} ', ' C8List{i,1}];
                    end
                end
            end
        end
    end
end
                
            

