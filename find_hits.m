function [ hit_list ] = find_hits( New_list,ref_list)
% Find links from New_list which match links in the ref_list, and store
% found links in output hit_list
% Yao Fu, May, 2010
count=0;
countn=0;
for i=1:size(New_list)
    hit=0;
    for j=1:size(ref_list)
        if strcmpi(New_list{i,1},ref_list{j,1}) && strcmpi(New_list{i,2},ref_list{j,2})
            hit=1;
            p=j;
        end
    end
    if hit==1;
         countn=countn+1;
         hit_list(countn,1:3)=New_list(i,:);
         hit_list(countn,4)=ref_list(p,4);
         if isempty(strfind(ref_list{p,4},New_list{i,3}))
            hit_list{countn,5}=0;
         else hit_list{countn,5}=1;
         end
    end
end
if countn==0;
    fprintf('%s\n','Sorry, no new regulation!');
    hit_list=[];
end    
end
    

