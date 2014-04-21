function [ z_M, cor_M ] = clr_cor(TFA,GeneE)
%	Borrelation based Context Likelihood Relatedness 
%   z_M: statistice likelihood score matrix
%   cor_M: Correlation score matrix
%   TFA: TFA matrix, GeneE: Gene expression matrix
%   by Fu, Yao(Al)
cor=corrcoef([TFA' GeneE']);
cor=cor.^2;
cor=cor-diag(diag(cor));
cor=cor(1:size(TFA),size(TFA)+1:size(cor));
cor_M=cor;
z_X=zscore(cor);
z_X(z_X<0)=0;
z_Y=zscore(cor');
z_Y(z_Y<0)=0;
z_M=sqrt(z_X.^2+z_Y'.^2);
end

