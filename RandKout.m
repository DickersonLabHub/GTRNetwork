function [ RC ] = RandKout( C,rate )
% Randomly knock out part of regulatory links from the connectivity matrix
%   C: connectivity matrix
%   rate: percenctage of links want to be knocked out, for 0 - 1
%  RC: connectivity matrix having some links randomly knocked out
% by Fu, Yao(Al)

num=int32(sum(sum(C))*(1-rate));
A=rand(size(C)).*C;
S=A(1,:);
for i=2:size(A,1)
    S=[S A(i,:)];
end
S=sort(S);
thr=S(size(S,2)-num);
RC=repmat(0,size(A));
RC(A>thr)=1;
end

