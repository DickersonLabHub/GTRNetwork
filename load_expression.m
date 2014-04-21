function [ Expression ] = load_expression( filename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
data=importdata(filename);
Expression.E=data.data;
Expression.Experiment_Name=data.textdata(1,2:end);
Expression.Gene_id=data.textdata(2:end,1);


end

