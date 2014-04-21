function [ summary ] = SumResults( GTRN_result )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[summary.APMI_CLR summary.APMI_CLRwOP]=SumTest(1,GTRN_result);
[summary.APMI_nonCLR summary.APMI_nonCLRwOP]=SumTest(2,GTRN_result);
[summary.correlation_CLR summary.correlation_CLRwOP]=SumTest(3,GTRN_result);    
[summary.correlation_nonCLR summary.correlation_nonCLRwOP]=SumTest(4,GTRN_result);    
end

function [result resultOP] = SumTest(row,GTRN_result)
for i=1:20
result(i,1)=size(GTRN_result(row,i).list,1);
result(i,2)=size(GTRN_result(row,i).NewLinks,1);
result(i,3)=size(GTRN_result(row,i).RefVerified,1);
if result(i,3)~=0
    result(i,4)=sum(cell2mat(GTRN_result(row,i).RefVerified(:,5)));
else
    result(i,4)=0;
end
result(i,5)=size(GTRN_result(row,i).BSVereified,1);
resultOP(i,1)=size(GTRN_result(row,i).Operon,1);
resultOP(i,2)=size(GTRN_result(row,i).NewLinksOP,1);
resultOP(i,3)=size(GTRN_result(row,i).RefVerifiedOP,1);
if resultOP(i,3)~=0
    resultOP(i,4)=sum(cell2mat(GTRN_result(row,i).RefVerifiedOP(:,5)));
else
    resultOP(i,4)=0;
end
resultOP(i,5)=size(GTRN_result(row,i).BSVereifiedOP,1);
end
end



