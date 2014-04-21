function [ new_exp_data ] = log10ratio( exp_data,base_exp_id )
%	Calculate log2 ratio of gene expression data based on control
%	conditions id 'base_exp_id'
%
%   by Fu, Yao(Al)
new_exp_data=exp_data;
if isnumeric(base_exp_id)
    b=base_exp_id;
else
    b=strmatch(base_exp_id,exp_data.Experiment_Name','exact');
end
if isempty(b)
    b=0;
    fprintf('%s\n','Error,no matching base experiment id');
else
    if b>size(exp_data.E,2)
        fprintf('%s\n','base_exp_id exceed number of experiments');
    else
    R=exp_data.E./repmat(exp_data.E(:,b),1,size(exp_data.E,2));
    new_exp_data.R=log10(R);
    end
end
end

