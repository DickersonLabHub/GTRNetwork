function [TFA] = SIMPLS( Data )
%	Partial Least Square to predict TFA
%   by Fu, Yao(Al)
E=Data.R;
A=Data.C;
[AL,EL,AS,ES,P,PCTVAR,MSE]=plsregress(A,E,rank(A),'cv',3,'mcreps',100);
[m,v]=min(MSE(2,:));
[A,EL,AS,ES,P]=plsregress(A,E,v-1);
TFA=P(2:size(P,1),:);
end

