function S = LoadRMA(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

data=importdata(filename);
S.E=data.data;
S.TimePoints=data.textdata(1,2:size(data.textdata,2));
S.GeneLabels=data.textdata(2:size(data.textdata,1),1);
end

