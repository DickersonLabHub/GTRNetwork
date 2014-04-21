function [ Prediction TFA Z M] = GTRNetwork( Expression, Regulatory,TFA_infer,Relatedness_infer,CLR,Network_size,organism,Operon,Binding_site)
%   GTRNetowk Algorithm to reconstruct gene regulatory networks
%   by Fu, Yao(Al)

% Predict TFA
TFA=GetTFA(Expression, Regulatory, TFA_infer,organism);
TFA=DirectTFA(TFA,Expression,Regulatory);

% Predict regulatory links
Old_list=list_net(Regulatory.C,0.5,Regulatory.TF_id,Regulatory.Gene_id);
if Network_size==0 % Testing model
    fprintf('%s\n','Network size is set to 0, testing model');
    fprintf('%s\n','Using relatedness score of Adaptive Partitioning Mutual Information ');
        [Z M]=clr_APMI(TFA.TFA,Expression.R);
    for i=1:20
        fprintf('%s%i\n','Reconstructing network at size ', i*100);
        thr=find_thr(i*100,Z);
        Prediction(1,i).list=list_net(Z,thr,TFA.tf_id,Expression.Gene_id); % 1: APMI with CLR
        Prediction(1,i).NewLinks=compare_lists(Prediction(1,i).list,Old_list);
        Prediction(1,i).NewLinks=DirectLinks(Prediction(1,i).NewLinks,TFA.TFA,TFA.tf_id(:,1),Expression.R,Expression.Gene_id(:,1));
        Prediction(1,i).RefVerified=find_hits(Prediction(1,i).NewLinks,Regulatory.AI);
        if exist('Binding_site')
            Prediction(1,i).BSVereified=check_bidingsite(Prediction(1,i).NewLinks,Binding_site,Operon);
        end
        if exist('Operon')
            Prediction(1,i).Operon=check_operon(Prediction(1,i).list,Operon);
            Prediction(1,i).NewLinksOP=compare_lists(Prediction(1,i).Operon,Old_list);
            Prediction(1,i).NewLinksOP=DirectLinks(Prediction(1,i).NewLinksOP,TFA.TFA,TFA.tf_id(:,1),Expression.R,Expression.Gene_id);
            Prediction(1,i).RefVerifiedOP=find_hits(Prediction(1,i).NewLinksOP,Regulatory.AI);
            if exist('Binding_site')
                Prediction(1,i).BSVereifiedOP=check_bidingsite(Prediction(1,i).NewLinksOP,Binding_site,Operon);
            end
        end
        thr=find_thr(i*100,M);
        Prediction(2,i).list=list_net(M,thr,TFA.tf_id,Expression.Gene_id); % 2: APMI w/o CLR
        Prediction(2,i).NewLinks=compare_lists(Prediction(2,i).list,Old_list);
        Prediction(2,i).NewLinks=DirectLinks(Prediction(2,i).NewLinks,TFA.TFA,TFA.tf_id(:,1),Expression.R,Expression.Gene_id);
        Prediction(2,i).RefVerified=find_hits(Prediction(2,i).NewLinks,Regulatory.AI);
        if exist('Binding_site')
            Prediction(2,i).BSVereified=check_bidingsite(Prediction(2,i).NewLinks,Binding_site,Operon);
        end
        if exist('Operon')
            Prediction(2,i).Operon=check_operon(Prediction(2,i).list,Operon);
            Prediction(2,i).NewLinksOP=compare_lists(Prediction(2,i).Operon,Old_list);
            Prediction(2,i).NewLinksOP=DirectLinks(Prediction(2,i).NewLinksOP,TFA.TFA,TFA.tf_id(:,1),Expression.R,Expression.Gene_id);
            Prediction(2,i).RefVerifiedOP=find_hits(Prediction(2,i).NewLinksOP,Regulatory.AI);
            if exist('Binding_site')
                Prediction(2,i).BSVereifiedOP=check_bidingsite(Prediction(2,i).NewLinksOP,Binding_site,Operon);
            end
        end
        fprintf('%s%i%s%i','APMI-CLR: NewLinks:',size(Prediction(1,i).NewLinks,1),' Verified:',size(Prediction(1,i).RefVerified,1));
        fprintf('%s%i%s%i\n',' NewLinksOP:',size(Prediction(1,i).NewLinksOP,1),' VerifiedOP:',size(Prediction(1,i).RefVerifiedOP,1));
        fprintf('%s%i%s%i','APMI: NewLinks:',size(Prediction(2,i).NewLinks,1),' Verified:',size(Prediction(2,i).RefVerified,1));
        fprintf('%s%i%s%i\n',' NewLinksOP:',size(Prediction(2,i).NewLinksOP,1),' VerifiedOP:',size(Prediction(2,i).RefVerifiedOP,1));
    end
    fprintf('%s\n','Using relatedness score of Pearson Correlation ');
        [Z M]=clr_cor(TFA.TFA,Expression.R);
    for i=1:20
        fprintf('%s%i\n','Reconstructing network at size ', i*100);
        thr=find_thr(i*100,Z);
        Prediction(3,i).list=list_net(Z,thr,TFA.tf_id,Expression.Gene_id); % 3: Correlation w/ CLR
        Prediction(3,i).NewLinks=compare_lists(Prediction(3,i).list,Old_list);
        Prediction(3,i).NewLinks=DirectLinks(Prediction(3,i).NewLinks,TFA.TFA,TFA.tf_id(:,1),Expression.R,Expression.Gene_id);
        Prediction(3,i).RefVerified=find_hits(Prediction(3,i).NewLinks,Regulatory.AI);
        if exist('Binding_site')
            Prediction(3,i).BSVereified=check_bidingsite(Prediction(3,i).NewLinks,Binding_site,Operon);
        end
        if exist('Operon')
            Prediction(3,i).Operon=check_operon(Prediction(3,i).list,Operon);
            Prediction(3,i).NewLinksOP=compare_lists(Prediction(3,i).Operon,Old_list);
            Prediction(3,i).NewLinksOP=DirectLinks(Prediction(3,i).NewLinksOP(:,1:2),TFA.TFA,TFA.tf_id(:,1),Expression.R,Expression.Gene_id);
            Prediction(3,i).RefVerifiedOP=find_hits(Prediction(3,i).NewLinksOP,Regulatory.AI);
            if exist('Binding_site')
                Prediction(3,i).BSVereifiedOP=check_bidingsite(Prediction(3,i).NewLinksOP,Binding_site,Operon);
            end
        end
        thr=find_thr(i*100,M);
        Prediction(4,i).list=list_net(M,thr,TFA.tf_id,Expression.Gene_id); % 4: Correlation w/o CLR
        Prediction(4,i).NewLinks=compare_lists(Prediction(4,i).list,Old_list);
        Prediction(4,i).NewLinks=DirectLinks(Prediction(4,i).NewLinks,TFA.TFA,TFA.tf_id(:,1),Expression.R,Expression.Gene_id);
        Prediction(4,i).RefVerified=find_hits(Prediction(4,i).NewLinks,Regulatory.AI);
        if exist('Binding_site')
            Prediction(4,i).BSVereified=check_bidingsite(Prediction(4,i).NewLinks,Binding_site,Operon);
        end
        if exist('Operon')
            Prediction(4,i).Operon=check_operon(Prediction(4,i).list,Operon);
            Prediction(4,i).NewLinksOP=compare_lists(Prediction(4,i).Operon,Old_list);
            Prediction(4,i).NewLinksOP=DirectLinks(Prediction(4,i).NewLinksOP,TFA.TFA,TFA.tf_id(:,1),Expression.R,Expression.Gene_id);
            Prediction(4,i).RefVerifiedOP=find_hits(Prediction(4,i).NewLinksOP,Regulatory.AI);
            if exist('Binding_site')
                Prediction(4,i).BSVereifiedOP=check_bidingsite(Prediction(4,i).NewLinksOP,Binding_site,Operon);
            end
        end
        fprintf('%s%i%s%i','Corr-CLR: NewLinks:',size(Prediction(3,i).NewLinks,1),' Verified:',size(Prediction(3,i).RefVerified,1));
        fprintf('%s%i%s%i\n',' NewLinksOP:',size(Prediction(3,i).NewLinksOP,1),' VerifiedOP:',size(Prediction(3,i).RefVerifiedOP,1));
        fprintf('%s%i%s%i','Correlation: NewLinks:',size(Prediction(4,i).NewLinks,1),' Verified:',size(Prediction(4,i).RefVerified,1));
        fprintf('%s%i%s%i\n',' NewLinksOP:',size(Prediction(4,i).NewLinksOP,1),' VerifiedOP:',size(Prediction(4,i).RefVerifiedOP,1));
    end
        
else %Prediction model
    if strcmpi(Relatedness_infer,'APMI')
        fprintf('%s\n','Using relatedness score of Adaptive Partitioning Mutual Information ');
        [Z M]=clr_APMI(TFA.TFA,Expression.R);
    elseif strcmpi(Relatedness_infer,'Correlation')
        fprintf('%s\n','Using relatedness score of Pearson Correlation ');
        [Z M]=clr_cor(TFA.TFA,Expression.R);
    else
        fprintf('%s\n','Please choose relatedness inference algorithm from APMI and Correlation');
        fprintf('%s\n','Here using correlation as default');
        [Z M]=clr_cor(TFA.TFA,Expression.R);
    end
    if CLR>0
        thr=find_thr(Network_size,Z);
        Prediction.list=list_net(Z,thr,TFA.tf_id,Expression.Gene_id);
    else
        thr=find_thr(Network_size,M);
        Prediction.list=list_net(M,thr,TFA.tf_id,Expression.Gene_id);
    end 
    Prediction.NewLinks=compare_lists(Prediction.list,Old_list);
    Prediction.NewLinks=DirectLinks(Prediction.NewLinks,TFA.TFA,TFA.tf_id(:,1),Expression.R,Expression.Gene_id);
    Prediction.RefVerified=find_hits(Prediction.NewLinks,Regulatory.AI);
    if exist('Binding_site')
        Prediction.BSVereified=check_bidingsite(Prediction.NewLinks,Binding_site,Operon);
    end
    if exist('Operon')
        Prediction.Operon=check_operon(Prediction.list,Operon);
        Prediction.NewLinksOP=compare_lists(Prediction.Operon,Old_list);
        Prediction.NewLinksOP=DirectLinks(Prediction.NewLinksOP,TFA.TFA,TFA.tf_id(:,1),Expression.R,Expression.Gene_id);
        Prediction.RefVerifiedOP=find_hits(Prediction.NewLinksOP,Regulatory.AI);
        if exist('Binding_site')
        Prediction.BSVereifiedOP=check_bidingsite(Prediction.NewLinksOP,Binding_site,Operon);
        end
    end
    fprintf('%s%i%s%i%','NewLinks:',size(Prediction.NewLinks,1),' Verified:',size(Prediction.RefVerified,1));
    fprintf('%s%i%s%i%\n',' NewLinksOP:',size(Prediction.NewLinksOP,1),' VerifiedOP:',size(Prediction.RefVerifiedOP,1));    
    Prediction.TFA=TFA;
    Prediction.Z=Z;
    Prediction.M=M;
end

end

