function SubData=Subnetwork(Data)


Ao=Data.C;
Eo=Data.R;
if isempty(Eo)
    E_input=0;
    SubData.R=[];
else
    E_input=1;
end
tfName=Data.tf_id;

[N,L]=size(Ao);
A=sign(Ao).*rand(N,L);
 
if nargin<2
    thres=2;
end

%Check condition 3. eliminate any gene regulated by more than M TFs
M=size(Data.R,2);
C=sum(Data.C');
keep_genes=find(C<M);
Eo=Eo(keep_genes,:);
A=A(keep_genes,:);
geneName=Data.gene_id(keep_genes,:);


%% 


%screening the network to eliminate TFs regulating less than 3 genes
%and coordinate TFs
[A,infoTF,infoGenes]=screenTF(A,thres);

[A,flag]=trimingA(A);
if flag==0
    Anet=A;
elseif flag==2
    fprintf('Sorry!  The network is not NCA compliant.');
    Anet=A;
    S=[];
else
    %make sure number of regulated genes is greater or equal to thres
    [A,infoTF,infoGenes]=rescreen(A,thres); 
    t1=find(infoTF);
    t2=find(infoGenes);
    At=A(t2,t1);
    h=CheckingRank(At);    
    while h==0 
        [A,temp_flag]=trimingA(A);
        if temp_flag==0
            warndlg('Sorry!  The network is not NCA compliant.');
            Anet=A;
            SubData=[];
            h=2;  %jump out of while loop
        else
            [A,infoTF,infoGenes]=rescreen(A,thres);
            t1=find(infoTF);
            t2=find(infoGenes);
            At=A(t2,t1);
            h=CheckingRank(At);
        end
    end
        
    if h==1
        Anet=A;
        col_sec=find(infoTF);
        row_sec=find(infoGenes);
        Asub=Anet(row_sec,col_sec);
        SubData.C=Asub;
        SubData.gene_id=geneName(row_sec,:);
        SubData.tf_id=tfName(col_sec,:);
        SubData.Mg=row_sec;
        SubData.Mtf=col_sec;
        if E_input==1
            SubData.R=Eo(row_sec,:);
            %S.exptname=exptName;
            %else
            %SubData.R=[];
            %S.exptname=[];
        end
        %d=printOutput(S,E_input);    
    end   
end
%SubData.tf_id=SubData.tf_id';
%SubData.gene_id=SubData.gene_id';

end
%----------------------------------------------
function [newA,tfInf,geneInf]=screenTF(oldA,thres)

[N,L]=size(oldA);
tfInf=ones(L,1);
geneInf=ones(N,1);
newA=oldA;
numGenes=sum(abs(sign(oldA)));  %row vector
[ord_num, tfIndices]=sort(numGenes);  %rows vectors

for k=1:L-1
    id=tfIndices(k);
    regg=find(newA(:,id)~=0);
    num_regg=length(regg);
    if num_regg < thres     %if number of regulated genes is less than 3
        newA(regg,:)=zeros(num_regg,L);
        tfInf(id)=0;
        geneInf(regg)=zeros(num_regg,1);
        fprintf('Note: The TF # %g regulates less than %g genes, so it was removed.\n',id,thres);
    elseif num_regg<=N-L+1
        temp=newA(regg,tfIndices(k+1:L));
        temp_prop=all(temp);   %which col(s) of temp has (have) ALL nonzero elements
        if any(temp_prop)       %if there exist(s) such column(s)
            tfInd(id)=0;
            geneInd(regg)=zeros(num_regg,1);
            newA(regg,:)=zeros(num_regg,L);
            parent_TF=k+find(temp_prop);
            fprintf('The TF # %g always co-regulates with: ',id);
            fprintf(' %g \n',tfIndices(parent_TF));
        end
    else
        %warndlg('There are too few genes per TF, so the network is not NCA compliant');
    end
        
end     
end
%-----------------------------------------
function [newA, flag]=trimingA(oldA)

L=size(oldA,2);
iter=1;
flag=0;
newA=oldA;
%Iteration loops modify network for NCA compliance
while (flag==0 & iter<1000)
    r_total=0;         %reset r_total
    tf_prop=any(newA);    %row vector
    chosenTF=find(tf_prop);      %row vector
    numTF=length(chosenTF);      %scalar
    
    if numTF==1  %don't need to check NCA if there is exist only 1 TF in subnetwork
        jump=1;
    else
        jump=0;
    end
    
    %Check NCA criteria
    k=1;
    while k<(numTF+1) & jump==0
        id=chosenTF(k);
        nregg=find(newA(:,id)==0);
        if k==1         %define the reduced matrix from the subnetwork of chosenTF
            tempTF=chosenTF(k+1:numTF);
            Ar=newA(nregg,tempTF);
        else
            tempTF=chosenTF([1:k-1 k+1:numTF]);
            Ar=newA(nregg,tempTF);
        end
        
        r=rank(Ar);        %Ar is nregg x numTF-1 
        r_total=r_total+r;
        
        if r<numTF-1      
            %Ar has at least 2 columns having only one nonzero element in the identical row
            %or Ar has zero column
            Ar=abs(sign(Ar));
            ng_perArcol=sum(Ar);
            bad_ind0=find(ng_perArcol==0);   %Ar has zero column
            if ~isempty(bad_ind0)
                badTF0=tempTF(bad_ind0);
                nbadTF0=length(badTF0);
                for jjj=1:nbadTF0
                    temp_id=badTF0(jjj); 
                    temp_genes=find(newA(:,temp_id));
                    temp_no=length(temp_genes);
                    newA(temp_genes,:)=zeros(temp_no,L);
                end
            end    
            bad_indAr=find(ng_perArcol==1);   %Ar has more than 1 column having only single nonzero element
            if ~isempty(bad_indAr)
                badTF=tempTF(bad_indAr); %vector stores id of 'bad' TFs
                nbad=length(badTF);
                ng_badTF=[];
                for j=1:nbad
                    temp_ng=length(find(newA(:,badTF(j))));
                    ng_badTF=[ng_badTF;temp_ng];
                end
                [temp_v,temp_order]=sort(ng_badTF);
                for jj=1:nbad-1
                    temp_id=badTF(temp_order(jj));
                    temp_genes=find(newA(:,temp_id));
                    temp_no=length(temp_genes);
                    newA(temp_genes,:)=zeros(temp_no,L);
                end
            end
            jump=1;   %exit the while loop checking NCA criteria
        end
        k=k+1;
    end
      
    if jump==1
        flag=0;
    else
        if r_total<numTF*(numTF-1)| r_total==0  %rank deficiency or numberTF=1
            flag=2;
        else
            flag=1;
        end
    end
    iter=iter+1;
end
end
%-----------------------------------------------
function [Apattern,finalTF,finalGenes]=rescreen(oldA,cutoff)

[N,L]=size(oldA);

finalTF=ones(L,1);
finalGenes=ones(N,1);
Apattern=sign(oldA);

nGenes=sum(abs(Apattern));

for k=1:L
    ng=nGenes(k);
    if ng<cutoff & ng~=0
        regg=find(oldA(:,k)~=0);
        Apattern(regg,:)=zeros(length(regg),L);
    end
end

nTf=sum(abs(Apattern),2);  %col vector
finalGenes=sign(nTf);
%nrGenes=find(nTf==0);
%finalGenes(nrGenes)=zeros(length(nrGenes),1);

nGenes=sum(abs(Apattern),1);  %row vector
finalTF=sign(nGenes);
%nrTFs=find(nGenes==0);
%finalTF(nrTFs)=zeros(length(nrTFs),1);  
end
%----------------------------------
function h = CheckingRank(alpha)

[N,L]= size(alpha);

alpha=rand(N,L).*sign(alpha);  % using random nonzero number instead of 1 and 0

R= 0;   % --> initialize rank of the augmented sparse system
nz= 0;  % --> initialize total # of zeros in the augmented system

for l= 1:L,
    v= find(alpha(:,l)==0);
    nz= nz+length(v);
    R= R+ rank(alpha(v,:));
end

%R
%pause
if (R < L*(L-1))
    h=0;
else
    h=1;
end
end

%----------------------------
function d=printOutput(S,inputE)

d=1;
D=SubData.C;
G=SubData.gene_id;
%G=char(G);
T=SubData.tf_id;
T=char(T);

if inputE==1
    E=SubData.R;
    Expt=S.exptname;
    Expt=char(Expt);
    M=size(E,2);
end
end

