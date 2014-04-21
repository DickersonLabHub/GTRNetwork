function [ link_direction ] = DirectLinks( link_list, TFA, TFA_id, Expression, gene_id )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
link_direction=link_list;
for i =1:size(link_list,1)
    if ~isempty(strmatch(link_list(i,2),gene_id,'exact'))
    if corrcoef(TFA(strmatch(link_list(i,1),TFA_id,'exact'),:),Expression(strmatch(link_list(i,2),gene_id,'exact'),:))>0
        link_direction{i,3}='+';
    else
        link_direction{i,3}='-';
    end
    else
        link_direction{i,3}='';
    end
    
end    
end

