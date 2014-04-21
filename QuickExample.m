%Estimate direction corrected TFAs
%1. Load Expression file which is formated for this tool
load([cd '/Data/E.coliC8.mat']); %Please format your expression file following this file being loaded
%2. Calculate Log2Ratio
Expression=E2R(ExpressionC8,1);
%3. Load Gene Regulatory Network Topology
TF_Gene=load_TF_gene([cd '/Data/RegulonDB_6.7/network_tf_gene.txt']);
%4. Caculate TFA
TFA = GetTFA(Expression, TF_Gene, 'SVD', 'E.coli'); 
%   Note:   use GetTFAY for Yeast
%           TFA infer method could be 'SVD'(fastest), 'EM'(best
%           performance) or 'PLS'(complete set of TFAs)
%5. Correct TFA direction
Directed_TFA = DirectTFA(TFA,Expression,TF_Gene);

%Reconstruct Gene Regulatory Network
%1. Load Expression file from database
Expression = load_m3d([cd '/Data/m3d_ecoli/avg_E_coli_v4_Build_6_exps466probes4297.tab'],'E.coli');
%   Note: use load_m3dY for Yeast
%2. Calculate Log2Ratio
Expression=E2R(Expression,110); 
%3. Load Gene Regulatory Network Topology
TF_Gene=load_TF_gene([cd '/Data/RegulonDB_6.7/network_tf_gene.txt']);
%   Note: use load_TF_gene_Yeast for Yeast date from eastract database
%Optional: load operon and binding site information
Operon=load_operon([cd '/Data/RegulonDB_6.7/OperonSet.txt']);
BindingSite=load_bindingsites([cd '/Data/RegulonDB_6.7/BindingSiteSet.txt']);
%4. Reconstruct
[Prediction TFA Z M] = GTRNetwork(Expression,TF_Gene,'SVD', 'correlation',1,100,'E.coli',Operon,BindingSite);
%   Note:   used GTRNetowrkY for yeast (NO need to specify organsim for Yeast using GTRNetworkY)
%           Prediction is the predicted regulatory links at the difined
%           network size£¨here is set at 100)
%           TFA is the predicted TFAs
%           Z are the likelihood score matrix
%           M is the caculated relatedness 
%           'SVD' can be replaced by 'EM' or 'PLS'
%           'correlation' can be replaced by 'APMI'(slower, but potential
%           have better performance)
%           Operon and binding_site are optional













