function [ list,count ] = list_net( M,thr,TFA_id,Gene_id)
%   List network links at certain threshold (thr) from a score matrix (M) 
%   by Fu, Yao£¨Al)
count=0;
for i=1:size(M,1)
    for j=1:size(M,2)
        if abs(M(i,j))>=thr
            count=count+1;
            list(count,:)=[TFA_id(i,1) Gene_id(j,1) M(i,j)];
        end
    end
end
if count==0
    list=[' ' ' '];
end

end

