function [ verif ] = check_bidingsite( link_list, bindingsite, operon )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n=1;
for i=1:size(operon.id,1)
    for j=1:size(operon.gene,2)
        if (operon.gene{i,j})
            Op(n,1)=operon.id(i,1);
            Op(n,2)=operon.gene(i,j);
            n=n+1;
        end
    end
end
n=1;
for i=1:size(link_list,1)
    t=strmatch(link_list(i,2),Op(:,2),'exact');
    if ~isempty(t)
        for j=1:size(t,1);
        p=strmatch(Op(t(j,1),1),bindingsite(:,2),'exact');
        sp=size(p,1);
        q=strmatch(link_list(i,1),bindingsite(p,1),'exact');
        if ~isempty(q)
            verif(n,:)=link_list(i,:);
            n=n+1;
        end
        end
    end
end
if ~exist('verif')
    fprintf('%s\n','No predicted links can be verified by TF binding sites.');
    verif='';
end
end

