function [TFA] = GetTFAY( Expression, Regulatory,TFA_infer)
%   GTRNetowk Algorithm to reconstruct gene regulatory networks
%   by Fu, Yao(Al)

if strcmpi(TFA_infer,'None')
    fprintf('%s\n','Using TF expression level as TFA ');
    N=1;
    for i=1:size(Regulatory.C,1)
        n=strmatch(Regulatory.TF_id(i,1),Expression.Gene_id(:,1),'exact');
        if ~isempty(n)
        TFA.TFA(N,:)=Expression.R(n,:);
        TFA.tf_id(N,:)=Regulatory.TF_id(i,:);
        TFA.C(N,:)=Regulatory.C(i,:);
        N=N+1;
        end
    end
    TFA.gene_id=Regulatory.Gene_id;
else
    Data0=MatchGenesY(Expression,Regulatory);
    if strcmpi(TFA_infer,'PLS')
        fprintf('%s\n','Using Partial Least Square as the TFA infer algorithm ');
        TFA.TFA=SIMPLS(Data0);
        Data1=Data0;
    else
        Data1=Subnetwork(Data0);
        if strcmpi(TFA_infer,'EM')
            fprintf('%s\n','Using generalized Network Component Analysis with regularizaion as the TFA infer algorithm ');
            [TFA CS]=gNCAr(Data1);
            [CS0 TFA.TFA]=Normalization(CS,TFA);
        elseif strcmpi(TFA_infer,'SVD')
            fprintf('%s\n','Using Fast Network Component Analysis as the TFA infer algorithm ');
            TFA.TFA=FastNCA(Data1);
        end
    end
    TFA.tf_id=Data1.tf_id;
    TFA.C=Data1.R/TFA.TFA;
    TFA.gene_id=Data1.gene_id;
end
if ~exist('TFA')
    fprintf('%s\n','Please choose TFA infer algorithm from None, PLS, EM and SVD');
    fprintf('%s\n','Here using TF expression level as TFA for default');
    N=1;
    for i=1:size(Regulatory.C,1)
        n=strmatch(Regulatory.TF_id(i,1),Expression.Gene_id(:,1),'exact');
        if ~isempty(n)
        TFA.TFA(N,:)=Expression.R(n,:);
        TFA.tf_id(N,:)=Regulatory.TF_id(i,:);
        TFA.C(N,:)=Regulatory.C(i,:);
        N=N+1;
        end
    end
    TFA.gene_id=Regulatory.Gene_id;
end


end

