function TFA_dist=GetDistribution(Expression,Regulatory,TFA_infer)
ExpressionR=Expression;
ExpressionR.R=shuffle(Expression.R);
TFAR=GetTFA(ExpressionR,Regulatory,TFA_infer);
TFA.R(:,1:size(TFAR.TFA,2))=TFAR.TFA;
for i=1:49
    fprintf('%i\t',i+1);
ExpressionR=Expression;
ExpressionR.R=shuffle(Expression.R);
TFAR=GetTFA(ExpressionR,Regulatory,TFA_infer);
TFA.R(:,size(TFA.R,2)+1:size(TFA.R,2)+size(TFAR.TFA,2))=TFAR.TFA;
end

TFA_dist=[mean(TFA.R,2) std(TFA.R,0,2)];
end
