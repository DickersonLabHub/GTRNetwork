function [ exp_data ] = load_m3d( filename,organism )
%	Read m3d gene expression file
%   by Fu, Yao(Al)
if strcmpi(organism,'E.coli')
raw=importdata(filename);
exp_data.E=raw.data;
exp_data.Experiment_Name=raw.textdata(1,2:end);
for i=2:size(raw.textdata,1)
    geneids=strread(raw.textdata{i,1},'%s','delimiter','_');
    exp_data.Gene_id(i-1,:)=geneids(1:2,1)';
end
elseif strcmpi(organism,'Yeast')
    exp_data=load_m3dY(filename);
else
    exp_data=[];
    fprintf('%s','Wrong Organism!');
end
end

