%match matrix
ntf=size(Regulatory.TF_id,1);
ng=size(Expression.Gene_id,1);
for i=1:ntf
    s=strmatch(Regulatory.TF_id(i,1),TFA.tf_id(:,1),'exact');
    if isempty(s)
        z_mat(i,:)=repmat(0,1,ng);
        m_mat(i,:)=repmat(0,1,ng);
    else
        z_mat(i,:)=Z(s,:);
        m_mat(i,:)=M(s,:);
    end
end
matched_genes=0;
for i=1:ng
    s=strmatch(Expression.Gene_id(i,1),Regulatory.Gene_id(:,1),'exact');
    if isempty(s)
        initial_C(:,i)=repmat(0,ntf,1);
        ref(:,i)=repmat(0,ntf,1);
        refC(:,i)=repmat({[]},ntf,1);
    else
        initial_C(:,i)=Regulatory.C(:,s);
        ref(:,i)=Regulatory.ref(:,s);
        refC(:,i)=Regulatory.refC(:,s);
        matched_genes=matched_genes+1;
    end
end
opm=repmat(0,size(Operon.n,1),ng);
for i=1:size(Operon.n,1)
    for j=1:Operon.n(i,1)
        s=strmatch(Operon.gene(i,j),Expression.Gene_id(:,1),'exact');
        if ~isempty(s)
            opm(i,s)=1;
        end
    end
end
op_mat=repmat(0,ng,ng);
for i=1:ng
    if sum(opm(:,i))~=0
        op_mat(i,:)=opm(opm(:,i)==1,:);
        op_mat(i,i)=0;
    end
end


%test performance
N_R=0;
n=0;
while(N_R<0.5)
n=n+1;
thr=find_thr(10*n,z_mat);
z_tmp=repmat(0,ntf,ng);
z_tmp(z_mat>=thr)=1;
z_op=z_tmp;
for i=1:ntf
    z_op(i,:)=z_tmp(i,:)+sum(op_mat(z_tmp(i,:)==1',:),1);
end
z_op(z_op>0)=1;

G_T=sum(sum(ref));
N_T=G_T-sum(sum(initial_C));
G_P=sum(sum(z_op));
N_P=numel(z_tmp(z_op-initial_C==1));
G_V=G_P-numel(z_tmp(z_op-ref==1));
N_V=N_P-numel(z_tmp(z_op-ref==1));
G_PC=G_V/G_P;
N_PC=N_V/N_P;
G_R=G_V/G_T;
N_R=N_V/N_T;
Performance(n,1:4)=[G_PC G_R N_PC N_R];
if n>=400
    break;
end
end

