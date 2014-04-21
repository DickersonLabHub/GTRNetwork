function [ Network ] = Prediction2Network( Old_Network, Prediction )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

g=size(Old_Network.Gene_id,1);
l=size(Old_Network.AI,1);
Network=Old_Network;
for i=1:size(Prediction.NewLinksOP,1)
    s=strmatch(Prediction.NewLinksOP(i,2),Network.Gene_id(:,1),'exact');
    if isempty(s)
        g=g+1;
        s=g;
        Network.Gene_id(s,1)=Prediction.NewLinksOP(i,2);
    end
    t=strmatch(Prediction.NewLinksOP(i,1),Network.TF_id(:,1),'exact');
    Network.AI(l+i,[1 2 4])=Prediction.NewLinksOP(i,:);
    Network.C(t,s)=1;
    Network.CD(t,s)=Prediction.NewLinksOP(i,3);
end
Network.ref=Network.C;
Network.refC=Network.CD;
end


