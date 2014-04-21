function [ GTRN_result_new ] = debug3( GTRN_result )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
GTRN_result_new=GTRN_result;
for i=1:size(GTRN_result,2)
    for j=1:size(GTRN_result(3,i).NewLinks,1)
        p=strmatch(GTRN_result(3,i).NewLinks(j,1),GTRN_result(3,i).NewLinksOP(:,1),'exact');
        q=strmatch(GTRN_result(3,i).NewLinks(j,2),GTRN_result(3,i).NewLinksOP(p,2),'exact');
        GTRN_result_new(3,i).NewLinks(j,3)=GTRN_result(3,i).NewLinksOP(q(1,1),3);
    end
    for j=1:size(GTRN_result(3,i).RefVerified,1)
        p=strmatch(GTRN_result(3,i).RefVerified(j,1),GTRN_result(3,i).RefVerifiedOP(:,1),'exact');
        q=strmatch(GTRN_result(3,i).RefVerified(j,2),GTRN_result(3,i).RefVerifiedOP(p,2),'exact');
        GTRN_result_new(3,i).RefVerified(j,[3 5])=GTRN_result(3,i).RefVerifiedOP(q(1,1),[3 5]);
    end
end

end

