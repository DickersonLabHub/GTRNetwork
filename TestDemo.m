%Test Demo
%Generate connectivity data by randomly knocking out links from latest
%RegulonDB links data
Expression_file=[cd '/Data/m3d_ecoli/avg_E_coli_v4_Build_6_exps466probes4297.tab'];
Regulatory_links=[cd '/Data/RegulonDB_6.7/network_tf_gene.txt'];
Operon_file=[cd '/Data/RegulonDB_6.7/OperonSet.txt'];
Binding_site_file=[cd '/Data/RegulonDB_6.7/BindingSiteSet.txt'];

Expression=load_m3d(Expression_file);
TF_Gene=load_TF_gene(Regulatory_links);
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
    g=0;
    while g==0
    TF_Gene50(i)=TF_Gene;
    TF_Gene50(i).C=RandKout(TF_Gene.C,0.5);
    g=checkgNCA(Expression,TF_Gene50(i));
    end 
    g=0;
    while g==0
    TF_Gene70(i)=TF_Gene;
    TF_Gene70(i).C=RandKout(TF_Gene.C,0.3);
    g=checkgNCA(Expression,TF_Gene70(i));
    end 
    g=0;
    while g==0
    TF_Gene90(i)=TF_Gene;
    TF_Gene90(i).C=RandKout(TF_Gene.C,0.1);
    g=checkgNCA(Expression,TF_Gene90(i));
    end
end
for i=1:5
    TF_Gene30(i).list=list_net(TF_Gene30(i).C,0.5,TF_Gene.TF_id,TF_Gene.Gene_id);
    TF_Gene50(i).list=list_net(TF_Gene50(i).C,0.5,TF_Gene.TF_id,TF_Gene.Gene_id);
    TF_Gene70(i).list=list_net(TF_Gene70(i).C,0.5,TF_Gene.TF_id,TF_Gene.Gene_id);
    TF_Gene90(i).list=list_net(TF_Gene90(i).C,0.5,TF_Gene.TF_id,TF_Gene.Gene_id);
end
%Run complete test
T30=Autotest(Expression,TF_Gene30,Operon,BindingSite);
T50=Autotest(Expression,TF_Gene50,Operon,BindingSite);
T70=Autotest(Expression,TF_Gene70,Operon,BindingSite);
T90=Autotest(Expression,TF_Gene90,Operon,BindingSite);