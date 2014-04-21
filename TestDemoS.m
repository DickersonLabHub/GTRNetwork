%Test Demo
%Generate connectivity data by randomly knocking out links from latest
%RegulonDB links data
Expression_file=[cd '/Data/m3d_ecoli/avg_E_coli_v4_Build_6_exps466probes4297.tab'];
Regulatory_links=[cd '/Data/RegulonDB_7.0/network_tf_gene.txt'];
Operon_file=[cd '/Data/RegulonDB_7.0/OperonSet.txt'];
Binding_site_file=[cd '/Data/RegulonDB_7.0/BindingSiteSet.txt'];
Sigma_file=[cd '/Data/RegulonDB_7.0/network_sigma_gene.txt'];

Expression=load_m3d(Expression_file);
TF_Gene=load_TF_gene(Regulatory_links,Sigma_file);
Operon=load_operon(Operon_file);
%Expression=log2ratio(Expression,110);
Expression=E2R(Expression,110);
BindingSite=load_bindingsites(Binding_site_file);

for i=1:5 % Generate Random Knock out data
    g=0;
    while g==0
    TF_Gene30(i)=TF_Gene;
    TF_Gene30(i).C=RandKout(TF_Gene.C,0.7);
    g=checkgNCA(Expression,TF_Gene30(i));
    end 
    TF_Gene30(i).CD=TF_Gene.CD;
    TF_Gene30(i).CD(TF_Gene30(i).C==0)={[]};
    g=0;
    while g==0
    TF_Gene50(i)=TF_Gene;
    TF_Gene50(i).C=RandKout(TF_Gene.C,0.5);
    g=checkgNCA(Expression,TF_Gene50(i));
    end 
    TF_Gene50(i).CD=TF_Gene.CD;
    TF_Gene50(i).CD(TF_Gene50(i).C==0)={[]};
    g=0;
    while g==0
    TF_Gene70(i)=TF_Gene;
    TF_Gene70(i).C=RandKout(TF_Gene.C,0.3);
    g=checkgNCA(Expression,TF_Gene70(i));
    end 
    TF_Gene70(i).CD=TF_Gene.CD;
    TF_Gene70(i).CD(TF_Gene70(i).C==0)={[]};
    g=0;
    while g==0
    TF_Gene90(i)=TF_Gene;
    TF_Gene90(i).C=RandKout(TF_Gene.C,0.1);
    g=checkgNCA(Expression,TF_Gene90(i));
    end
    TF_Gene90(i).CD=TF_Gene.CD;
    TF_Gene90(i).CD(TF_Gene90(i).C==0)={[]};
end
for i=1:5
    TF_Gene30(i).list=list_net(TF_Gene30(i).C,0.5,TF_Gene.TF_id,TF_Gene.Gene_id);
    TF_Gene50(i).list=list_net(TF_Gene50(i).C,0.5,TF_Gene.TF_id,TF_Gene.Gene_id);
    TF_Gene70(i).list=list_net(TF_Gene70(i).C,0.5,TF_Gene.TF_id,TF_Gene.Gene_id);
    TF_Gene90(i).list=list_net(TF_Gene90(i).C,0.5,TF_Gene.TF_id,TF_Gene.Gene_id);
end
%Run complete test
mkdir result;
AutotestS(Expression,TF_Gene30,Operon,BindingSite,'result/RS30');
AutotestS(Expression,TF_Gene50,Operon,BindingSite,'result/RS50');
AutotestS(Expression,TF_Gene70,Operon,BindingSite,'result/RS70');
AutotestS(Expression,TF_Gene90,Operon,BindingSite,'result/Rs90');