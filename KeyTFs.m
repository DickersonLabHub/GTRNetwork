function [ TFs ] = KeyTFs( Expression1, TFA1,Expression2,TFA2,Network,Gene_id,TF_id,Gene_p_value,TFA_p_value )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Gene_diff=mean(Expression1,2)-mean(Expression2,2);
Gene_diff(Gene_diff<0)=-1;
Gene_diff(Gene_diff>=0)=1;
Gene_diff(mattest(Expression1,Expression2)>Gene_p_value)=0;
n=1;
m=1;
for i=1:size(Gene_diff,1)
    if Gene_diff(i,1)==1
        UPGene(n,:)=Gene_id(i,:);
        n=n+1;
    elseif Gene_diff(i,1)==-1
        DownGene(m,:)=Gene_id(i,:);
        m=m+1;
    end
end
TFA_diff=mean(TFA1,2)-mean(TFA2,2);
TFA_diff(TFA_diff<0)=-1;
TFA_diff(TFA_diff>=0)=1;
TFA_diff(mattest(TFA1,TFA2)>TFA_p_value)=0;
n=1;
m=1;
for i=1:size(TFA_diff,1)
    if TFA_diff(i,1)==1
        TFs.Activated_TF(n,:)=TF_id(i,:);
        n=n+1;
    elseif TFA_diff(i,1)==-1
        TFs.Silented_TF(m,:)=TF_id(i,:);
        m=m+1;
    end
end

n=1;
 for i=1:size(TFs.Activated_TF,1)
     s=strmatch(TFs.Activated_TF(i,1),Network(:,1),'exact');
     if ~isempty(s)
        m=0;
        for j=1:size(s,1)
            if strmatch('+',Network(s(j,1),3),'exact')
               if ~isempty(strmatch(Network(s(j,1),2),UPGene,'exact'))
                   m=m+1;
               end
            elseif strmatch('-',Network(s(j,1),3),'exact')
                if ~isempty(strmatch(Network(s(j,1),2),DownGene,'exact'))
                    m=m+1;
                end
            else
                if ~isempty(strmatch(Network(s(j,1),2),[UPGene;DownGene],'exact'))
                    m=m+1;
                end
            end
        end
        if m/size(s,1)>0.5
            TFs.Key_TF(n,1)=TFs.Activated_TF(i,1);
            n=n+1;
        end
     end
 end
 for i=1:size(TFs.Silented_TF,1)
     s=strmatch(TFs.Silented_TF(i,1),Network(:,1),'exact');
     if ~isempty(s)
        m=0;
        for j=1:size(s,1)
            if strmatch('-',Network(s(j,1),3),'exact')
               if ~isempty(strmatch(Network(s(j,1),2),UPGene,'exact'))
                   m=m+1;
               end
            elseif strmatch('+',Network(s(j,1),3),'exact')
                if ~isempty(strmatch(Network(s(j,1),2),DownGene,'exact'))
                    m=m+1;
                end
            else
                if ~isempty(strmatch(Network(s(j,1),2),[UPGene;DownGene],'exact'))
                    m=m+1;
                end
            end
        end
        if m/size(s,1)>0.5
            TFs.Key_TF(n,1)=TFs.Silented_TF(i,1);
            n=n+1;
        end
     end
 end
TFs.UpGenes=UPGene;
TFs.DownGenes=DownGene;
end

