function [Test]=Autotest( Expression,Regulatory,Operon,BindingSite)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i=1:size(Regulatory,2)
    fprintf('%s%i\n','Testing round',i);
    Test(i).PnoneCC=GTRNetwork(Expression,Regulatory(i),'None','Correlation',1,1,Operon,BindingSite);
    Test(i).PsvdCC=GTRNetwork(Expression,Regulatory(i),'SVD','Correlation',1,1,Operon,BindingSite);
    Test(i).PemCC=GTRNetwork(Expression,Regulatory(i),'EM','Correlation',1,1,Operon,BindingSite);
    Test(i).PplsCC=GTRNetwork(Expression,Regulatory(i),'pls','Correlation',1,1,Operon,BindingSite);
    
    Test(i).PnoneCN=GTRNetwork(Expression,Regulatory(i),'None','Correlation',0,1,Operon,BindingSite);
    Test(i).PsvdCN=GTRNetwork(Expression,Regulatory(i),'SVD','Correlation',0,1,Operon,BindingSite);
    Test(i).PemCN=GTRNetwork(Expression,Regulatory(i),'EM','Correlation',0,1,Operon,BindingSite);
    Test(i).PplsCN=GTRNetwork(Expression,Regulatory(i),'pls','Correlation',0,1,Operon,BindingSite);
    
    Test(i).PnoneAC=GTRNetwork(Expression,Regulatory(i),'None','APMI',1,1,Operon,BindingSite);
    Test(i).PsvdAC=GTRNetwork(Expression,Regulatory(i),'SVD','APMI',1,1,Operon,BindingSite);
    Test(i).PemAC=GTRNetwork(Expression,Regulatory(i),'EM','APMI',1,1,Operon,BindingSite);
    Test(i).PplsAC=GTRNetwork(Expression,Regulatory(i),'pls','APMI',1,1,Operon,BindingSite);
    
    Test(i).PnoneAN=GTRNetwork(Expression,Regulatory(i),'None','APMI',0,1,Operon,BindingSite);
    Test(i).PsvdAN=GTRNetwork(Expression,Regulatory(i),'SVD','APMI',0,1,Operon,BindingSite);
    Test(i).PemAN=GTRNetwork(Expression,Regulatory(i),'EM','APMI',0,1,Operon,BindingSite);
    Test(i).PplsAN =GTRNetwork(Expression,Regulatory(i),'pls','APMI',0,1,Operon,BindingSite);
end
end

