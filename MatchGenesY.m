function [ Matched_data ] = MatchGenesY( exp,TF_gene )
%	Match genes of gene expression data and TF-gene network topology data
%   by Fu, Yao(Al)

[nt ngt]=size(TF_gene.C);
[ng ne]=size(exp.R);
C0=zeros(ng,nt);
R0=zeros(ng,ne);

for i=1:ngt
    n=strmatch(upper(TF_gene.Gene_id(i,1)),exp.Gene_id(:,1),'exact');
    if ~isempty(n)
        C0(n(1,1),:)=TF_gene.C(:,i)';
        R0(n,:)=exp.R(n,:);
    end
end

Mtf=(1:nt)'.*(sum(C0)'>0);
Matched_data.Mtf=Mtf(Mtf~=0);
Mg=(1:ng)'.*(sum(C0,2)>0);
Matched_data.Mg=Mg(Mg~=0);
Matched_data.tf_id=TF_gene.TF_id(Mtf>0,:);
Matched_data.gene_id=exp.Gene_id(Mg>0,:);
Matched_data.R=R0(sum(C0,2)>0,:);
C0=C0(sum(C0,2)>0,:);
Matched_data.C=C0(:,sum(C0,1)'>0);
end

