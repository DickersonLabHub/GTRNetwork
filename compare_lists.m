function [ new_list ] = compare_lists( New_list,old_list)
% Compare two list of regulatory links
%   find new links from New_list which is not in old_list and store new
%   links in output new_list
% Yao Fu, May, 2010
count=0;
countn=0;
for i=1:size(New_list)
    hit=0;
    for j=1:size(old_list)
        if strcmpi(New_list{i,1},old_list{j,1}) && strcmpi(New_list{i,2},old_list{j,2})
            hit=1;
        end
    end
    if hit==0;
         countn=countn+1;
         new_list(countn,:)=New_list(i,:);
    end
end
if countn==0;
    fprintf('%s\n','Sorry, no new regulation!');
    new_list=[];
end    
end
    

