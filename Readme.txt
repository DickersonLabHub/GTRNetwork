GTRNetwork Package modified for CBiRC use
Some useful functions:

GetTFA( Expression, Regulatory,TFA_infer,organism)
Returns TFA prediction of E.coli, organism option please choose "E.coli", TFA_infer can be "EM", "SVD", or "PLS"
Expreesion is the data loaded from M3D database file using load_m3d function, for other expression data, please follow the structure of loaded m3d data. 
E.coliC8.mat is the formated Expression data for CBiRC E.coli C8 stress experiments, could be loaded directly
Regulatory is the data loaded from RegulonDB database file using load_TF_gene function, for other regulatory toplogies, please follow the format of loaded Regulatory data.

GetTFAY is the same TFA prediction function specific to Yeast data

GTRNetwork is the functions to reconstruct regulatory networks
Returns predicted new regulatory links, predicted TFAs, CLR processed relatedness score matrix and the orignial relatedness score (correlation or mutual information) matrix
Expression,Regulatory,and TFA_infer options are the same as those for the funcion GetTFA.
Relatedness_infer options could be "APMI" or "Correlation"
CLR options could be 0 or 1 (don't use or use CLR)
Networ_size is defind for the exptected predicted network size (including recalled known links)
operon and binding_site are used to specify additional operon and binding_site information file to help predictions. (RegulonDB format)

GTRNetworkY is the function specific to Yeast.

More examples please see TestDemo and Autotest files