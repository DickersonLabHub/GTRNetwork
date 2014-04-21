function [z_M,MI ] = clr_APMI(TFA,GeneE)
%	APMI based Context Likelihood Relatedness 
%   z_M: statistice likelihood score matrix
%   MI: Mutual information score matrix
%   TFA: TFA matrix, GeneE: Gene expression matrix
%   by Fu, Yao(Al)
sizeT=size(TFA,1);
sizeG=size(GeneE,1);
for i=1:sizeT
    for j=1:sizeG
        MI(i,j)=APMI(TFA(i,:),GeneE(j,:));
    end
    fprintf('%i%s%i\n',i,'of', sizeT);
end
z_X=zscore(MI);
z_X(z_X<0)=0;
z_Y=zscore(MI');
z_Y(z_Y<0)=0;
z_M=sqrt(z_X.^2+z_Y'.^2);
sizeZ1=size(z_M,1);
sizeZ2=size(z_M,2);
end