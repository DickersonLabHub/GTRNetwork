function [ MI ] = APMI( X, Y )
%Adaptive Partitioning Mutual Information
%  Input: X,Y the two matrices to calculate the Mutual information. 
%   Output: MI the mutual information matrix
%   Algorithm developed by Liang K-C et al. 2008 (Gene Regulatory Network
%   Reconstruction Using Conditional Mutual Information. In. Edited by Wang
%   X, vol. 2008. EURASIP Journal on Bioinformatics and Systems Biology; 2008: 14 pages.)
% Implimented by Yao Fu for GTRNetwork algorithm, May, 2010
RX=sort(X);
RY=sort(Y);
CHI=chi2inv(0.95,3);
H=[min(X) max(X) min(Y) max(Y)];
P=[1 1 1];
total_n=size(X,2);
s=1;
s_old=0;
SX=repmat(0,1,total_n);
SY=repmat(0,1,total_n);
while (s~=s_old)
    s_old=s;
    A=H;
    B=P;
    p=1;
    for i=1:s
        a=median(RX(sum(X<A(i,1))+1:sum(X<A(i,2))));
        b=median(RY(sum(Y<A(i,3))+1:sum(Y<A(i,4))));
        SX(X<a)=-1;
        SX(X>a)=1;
        SX(X<A(i,1))=0;
        SX(X>A(i,2))=0;
        SY(Y<b)=-1;
        SY(Y>b)=1;
        SY(Y<A(i,3))=0;
        SY(Y>A(i,4))=0;
        N(1)=sum(SX+SY==-2);
        N(2)=sum(SX-SY==2);
        N(3)=sum(SX-SY==-2);
        N(4)=sum(SX+SY==2);
        N(5)=sum(N(1:4));
        if CHI< sum((N(5)/4-N(1:4)).^2/(N(5)/4))
            BB=[N(1) sum(SX==-1) sum(SY==-1); N(2) sum(SX==1) sum(SY==-1); N(3) sum(SX==-1) sum(SY==1); N(4) sum(SX==1) sum(SY==1)]/total_n;
            AA=[A(i,1) a A(i,3) b; a A(i,2) A(i,3) b; A(i,1) a b A(i,4); a A(i,2) b A(i,4)];
            H=[H(1:p-1,:);AA;H(p+1:s,:)];
            P=[P(1:p-1,:);BB;P(p+1:s,:)];
            p=p+4;
            s=s+3;
        else p=p+1;
        end        
    end
end
L=log2(P(:,1)./(P(:,2).*P(:,3)));
L(P(:,1)==0)=0;
MI=sum(P(:,1).*L);

end
