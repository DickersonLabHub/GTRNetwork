function [ thr ] = find_thr( num,A )
%	find threshold of a matrix A to so that certain number (num) of
%	elements are larger than the threshold
%   by Fu, Yao(Al)
S=A(1,:);
for i=2:size(A,1)
    S=[S A(i,:)];
end
S=abs(S);
S=sort(S);
thr=S(size(S,2)-num+1);
end

