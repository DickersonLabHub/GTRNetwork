% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
% This code is for normalization of CS(a) and TFA(p) matrix after NCA run
% [s,a,p]=NCA(E,A,0);
% L : the number of transcription factor
% N : the number of genes
% Usage: function [nor_a, nor_p]=normalization(a,p)
% nor_a: normalized a(CS); nor_p: normalized p(TFA)
% written by Young Lyeol Yang


function [nor_a, nor_p, nor_pdev]=Normalization(a,p,p_dev)
[N,L]=size(a);
res=[];
NF=[];
for(k=1:1:L)
    for(kk=1:1:N)
        if (a(kk,k)~=0)
            temp=abs(a(kk,k));
            res=[res;temp];
        end
    end
    temp1=mean(res,1);
    NF=[NF;temp1];
    res=[];
end

nor_pdev=zeros(size(p));

% The following codes are for calculating the normalization factor of CS columns
NF1=NF.^-1;

% The following codes are for calculating the normalized CS(a) and TFA(p)
for(k=1:1:L)
   nor_a(:,k)=a(:,k)*NF1(k,1);
   nor_p(k,:)=p(k,:)*NF(k,1);
   if (nargin > 2),
       nor_pdev(k,:)=p_dev(k,:)*NF(k,1);
   end
end
   