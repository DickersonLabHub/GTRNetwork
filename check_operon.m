function [ operon_list ] = check_operon( link_list, operon )
%   find genes in the same operon as the target genes and generate links
%   between these same operon genes and their regulators
%   by Fu, Yao(Al)
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
n_tf=1;
TF(n_tf,1)=link_list(1,1);
for i=2:size(link_list,1)
    if isempty(strmatch(link_list(i,1),TF,'exact'))
        n_tf=n_tf+1;
        TF(n_tf,1)=link_list(i,1);
    end
end
n_oplist=0;
for i=1:n_tf
    list=strmatch(TF(i,1),link_list(:,1),'exact');
    gene=link_list(list,2);
    o=[];
    for j=1:size(gene,1)
        p=strmatch(gene(j,1),Op(:,2),'exact');
        if isempty(p)
            n_oplist=n_oplist+1;
            operon_list(n_oplist,1)=TF(i,1);
            operon_list(n_oplist,2)=gene(j,1);
        else
        o=[o;Op(p,1)];
        end
    end
    if ~isempty(o)
    n_op=1;
    op_list(n_op,1)=o(1,1);
    for j=1:size(o,1)
        if isempty(strmatch(o(j,1),op_list(1:n_op,1),'exact'))
                n_op=n_op+1;
                op_list(n_op,1)=o(j,1);
        end
    end
    for j=1:n_op
        p=strmatch(op_list(j,1),Op(:,1),'exact');
        np=size(p,1);
        operon_list(n_oplist+1:n_oplist+np,1)=repmat(TF(i,1),np,1);
        operon_list(n_oplist+1:n_oplist+np,2)=Op(p,2);
        n_oplist=n_oplist+np;
    end
    end
        
end


end

