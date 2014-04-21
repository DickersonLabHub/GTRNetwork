function [ PerformSum ] = AssessPerformance(Expression,Regulatory,Result,Operon)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for i=1:size(Result,2)
    [PerformSum(i).PnoneCC.Perofmance PerformSum(i).PnoneCC.AUC]=testperform(Expression,Result(1,i).PnoneCC,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PsvdCC.Perofmance PerformSum(i).PsvdCC.AUC]=testperform(Expression,Result(1,i).PsvdCC,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PemCC.Perofmance PerformSum(i).PemCC.AUC]=testperform(Expression,Result(1,i).PemCC,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PplsCC.Perofmance PerformSum(i).PplsCC.AUC]=testperform(Expression,Result(1,i).PplsCC,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PnoneCN.Perofmance PerformSum(i).PnoneCN.AUC]=testperform(Expression,Result(1,i).PnoneCN,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PsvdCN.Perofmance PerformSum(i).PsvdCN.AUC]=testperform(Expression,Result(1,i).PsvdCN,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PemCN.Perofmance PerformSum(i).PemCN.AUC]=testperform(Expression,Result(1,i).PemCN,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PplsCN.Perofmance PerformSum(i).PplsCN.AUC]=testperform(Expression,Result(1,i).PplsCN,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PnoneAC.Perofmance PerformSum(i).PnoneAC.AUC]=testperform(Expression,Result(1,i).PnoneAC,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PsvdAC.Perofmance PerformSum(i).PsvdAC.AUC]=testperform(Expression,Result(1,i).PsvdAC,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PemAC.Perofmance PerformSum(i).PemAC.AUC]=testperform(Expression,Result(1,i).PemAC,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PplsAC.Perofmance PerformSum(i).PplsAC.AUC]=testperform(Expression,Result(1,i).PplsAC,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PnoneAN.Perofmance PerformSum(i).PnoneAN.AUC]=testperform(Expression,Result(1,i).PnoneAN,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PsvdAN.Perofmance PerformSum(i).PsvdAN.AUC]=testperform(Expression,Result(1,i).PsvdAN,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PemAN.Perofmance PerformSum(i).PemAN.AUC]=testperform(Expression,Result(1,i).PemAN,Regulatory(1,i),Operon,0.15,0.15);
    [PerformSum(i).PplsAN.Perofmance PerformSum(i).PplsAN.AUC]=testperform(Expression,Result(1,i).PplsAN,Regulatory(1,i),Operon,0.15,0.15);
    

end

