function AutotestS( Expression,Regulatory,Operon,BindingSite,filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i=1:size(Regulatory,2)
    fprintf('%s%i\n','Testing round',i);
    Psvd=GTRNetwork(Expression,Regulatory(i),'SVD','Correlation',1,0,Operon,BindingSite);
    fprintf('%s%i%s%s%s%i%s\n','Saving Psvd ',i,' to ',filename,'_svd',i,'.mat');
    save([filename '_svd' num2str(i) '.mat'],'Psvd');
    Pem=GTRNetwork(Expression,Regulatory(i),'EM','Correlation',1,0,Operon,BindingSite);
    fprintf('%s%i%s%s%s%i%s\n','Saving Pem ',i,' to ',filename,'_em',i,'.mat');
    save([filename '_em' num2str(i) '.mat'],'Pem');
    Ppls=GTRNetwork(Expression,Regulatory(i),'pls','Correlation',1,0,Operon,BindingSite);
    fprintf('%s%i%s%s%s%i%s\n','Saving Ppls ',i,' to ',filename,'_pls',i,'.mat');
    save([filename '_pls' num2str(i) '.mat'],'Ppls');
end

end

