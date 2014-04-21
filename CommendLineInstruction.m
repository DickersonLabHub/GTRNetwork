% Step 1: Choose organism: 'E.coli' or 'Yeast'
organism='E.coli';
% Or
organism='Yeast';

% Step 2: load transcriptome data
% load M3D data
% E.coli
filename=[cd '/Data/m3d_ecoli/avg_E_coli_v4_Build_6_exps466probes4297.tab'];
% Yeast
filename=[cd '/Data/Yeast/S_cerevisiae_v3_Build_1_norm/Yeast_Genome_S98_v2_norm/avg_Yeast_Genome_S98_v2_exps247probes5520.tab'];
Expression=load_m3d(filename,organism);
% Or load your own data
%Example of E.coli C8 fatty acid stress data (Please make your own data
%follow the same structure as the example matlab data!)
filename=[cd '/Data/E.coliC8.mat'];
Expression=load(filename);
%Or
filename=[cd '/Data/E.coliC8.txt'];
Expression=load_expression(filename);
%log2 ratio conversion
% If the transcriptome data is already in log 2 ratio
    Expression.R=Expression.E;
% If the transcriptome data is no in log 2 ratio
    c=1; %the column # of one control experiment
    Expression=log2ratio(Expression,c);

% Step 3: load initial regulatory data
% E.coli
    %RegulonDB
    filename=[cd '/Data/RegulonDB_7.0/network_tf_gene.txt'];
    Regulatory=load_TF_gene(filename);
    %Or RegulonDB + Predicted Links
        % 125 new links
    filename=[cd '/Data/Ecoli125Links.mat'];
        % 381 new links
    filename=[cd '/Data/Ecoli381Links.mat'];
    load(filename);
    
% Yeast
    %YeastTract
    filename=[cd '/Data/Yeast/YeastRegulation20101213.tsv'];
    Regulatory=load_TF_Gene_Yeast(filename);
    %Or YeastTract = Predicted Links
        % new links at detected network size at 100
    filename=[cd '/Data/YeastLinks100.mat'];
        % new links at detected network size at 200
    filename=[cd '/Data/YeastLInks200.mat'];
    load(filename);
    
% Step 4A: Predict TFAs
    %Choose TFA predition algoritm
    %'PLS'=partial least square, 'SVD'=FastNCA, 'EM'=gNca-r,
    %'None'=gene expression of TF genes (E.coli only)
    TFAmethod='SVD';
    %Run prediction
    TFA=GetTFA(Expression,Regulatory,TFAmethod,organism);
    
% Step 4A-2: Analysis TFA results
    %Key TFs in response to condition change
    %Take the example of E.coli C8 stress data, column 1-2
    %are experiment condition and column 3-5 are control condition
    %Define p-value threshold
    Gene_p=0.05;
    TFA_p=0.05;
    KeyTF=KeyTFs(Expression.E(:,1:2),TFA.TFA(:,1:2),Expression.E(:,3:5),TFA.TFA(:,3:5),Regulatory.AI(:,[1 2 4]),Expression.Gene_id,TFA.tf_id,Gene_p,TFA_p);
    %Find regulatory network of condition change
    %Map to pathway (E.coli Only)
    % See note.m in pathwayOfGenes folder
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4B: Reconstruct gene regulatory networks
    %Choose TFA predition algorithm
    %'PLS'=partial least square, 'SVD'=FastNCA, 'EM'=gNca-r,
    %'None'=gene expression of TF genes (E.coli only)
    TFAmethod='SVD';
    %Choose network reconstruction algorithm
    Relevance='APMI';
    Relevance='Correlation';
    %Choose network size you want to reconstruct
    Network_size=100;
    %Run GTRNetwork
        %E.coli
        %load Operon and Bindingsite information for E. coli
        Operon_file=[cd '/Data/RegulonDB_7.0/OperonSet.txt'];
        Binding_site_file=[cd '/Data/RegulonDB_7.0/BindingSiteSet.txt'];
        Operon=load_operon(Operon_file);
        BindingSite=load_bindingsites(Binding_site_file);
        %Run
        [ Prediction TFA Z M] = GTRNetwork( Expression, Regulatory,TFAmethod,Relevance,1,Network_size,Operon,BindingSite);
        %Yeast
        [ Prediction TFA Z M] = GTRNetworkY( Expression, Regulatory,TFAmethod,Relevance,1,Network_size);
    
    
% Step 4B-2: Export new regulatory links
    % Put new links into initial regulatory network
    Regulatory=Prediction2Network(Regulatory,Prediction);
    



