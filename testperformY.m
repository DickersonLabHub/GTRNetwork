function [Performance] = testperformY(Expression,Prediction,Regulatory)

Z=Prediction.Z;
%M=Prediction.M;
TFA=Prediction.TFA;

%match matrix
ntf=size(Regulatory.TF_id,1);
ng=size(Expression.Gene_id,1);
for i=1:ntf
    s=strmatch(Regulatory.TF_id(i,1),TFA.tf_id(:,1),'exact');
    if isempty(s)
        z_mat(i,:)=repmat(0,1,ng);
%        m_mat(i,:)=repmat(0,1,ng);
    else
        z_mat(i,:)=Z(s,:);
%        m_mat(i,:)=M(s,:);
    end
end
matched_genes=0;
for i=1:ng
    s=strmatch(Expression.Gene_id(i,1),upper(Regulatory.Gene_id(:,1)),'exact');
    if isempty(s)
        initial_C(:,i)=repmat(0,ntf,1);
        ref(:,i)=repmat(0,ntf,1);
    else
        initial_C(:,i)=Regulatory.C(:,s);
        ref(:,i)=Regulatory.ref(:,s);
        matched_genes=matched_genes+1;
    end
end
%test performance
for n=1:100
thr=find_thr(10*n,z_mat);
z_tmp=repmat(0,ntf,ng);
z_tmp(z_mat>=thr)=1;
%z_op=z_tmp;

%G_T=sum(sum(ref)); % G for Global, N for New
%N_T=G_T-sum(sum(initial_C));
%G_P=sum(sum(z_op));
%N_P=numel(z_tmp(z_op-initial_C==1));
%G_V=G_P-numel(z_tmp(z_op-ref==1));
%N_V=N_P-numel(z_tmp(z_op-ref==1));
%G_PC=G_V/G_P;
%N_PC=N_V/N_P;
%G_R=G_V/G_T;
%N_R=N_V/N_T;
Performance(n,1)=numel(z_tmp((z_tmp-initial_C==1)));

%AUC
end


end
