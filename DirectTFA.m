function [ Directed_TFA ] = DirectTFA( TFA, Expression, Regulatory )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
m=size(TFA.TFA,1);
l=size(Regulatory.CD,2);
D=repmat(0,m,l);
for i=1:m
    t=strmatch(TFA.tf_id(i,1),Regulatory.TF_id(:,1),'exact');
    for j=1:l
        if strcmp(Regulatory.CD(t,j),'-')
            D(i,j)=-1;
        elseif strcmp(Regulatory.CD(t,j),'+')
            D(i,j)=1;
        end
    end
end
for i=1:m
    CD=0;
    for j=1:l
        if D(i,j)~=0
            t=strmatch(Regulatory.Gene_id(j,1),Expression.Gene_id(:,1),'exact');
            if ~isempty(t)
            c=corrcoef(Expression.R(t,:),TFA.TFA(i,:));
            CD=CD+c(1,2)/D(i,j);
            end
        end
    end
    if CD<0
        Directed_TFA.TFA(i,:)=-TFA.TFA(i,:);
        Directed_TFA.C(:,i)=-TFA.C(:,i);
        Directed_TFA.CD(i)=CD;
    else Directed_TFA.TFA(i,:)=TFA.TFA(i,:);
        Directed_TFA.C(:,i)=TFA.C(:,i);
         Directed_TFA.CD(i)=CD;
    end
end
Directed_TFA.tf_id=TFA.tf_id;
Directed_TFA.D=D;

end

