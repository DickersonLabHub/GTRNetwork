function matrix= loadTomatrix(prob,exp,col,matrx) 
matrix=matrx;
for i=1:size(prob,1)
if ~isempty(exp(exp(:,1)==strread(prob{i,1},'%f',1),2))
matrix(i,col)=exp(exp(:,1)==strread(prob{i,1},'%f',1),2);
end
end