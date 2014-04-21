for i=1:6
for j=1:size(list(i,1).list,1)
if isempty(strmatch(list(i,1).list(j,2),TF_Gene.AI(:,2),'exact'))
fprintf(fid,'%s\t%s%s\t%i\t',list(i,1).list{j,1},list(i,1).list{j,2},'*',i*100);
else
fprintf(fid,'%s\t%s\t%i\t',list(i,1).list{j,1},list(i,1).list{j,2},i*100);
end
s=strmatch(list(i,1).list(j,2),EcoProb,'exact');
if ~isempty(s)
for k=1:size(Eco(1,s).str,1)
end
for k=2:size(Eco(1,s).str,1)
fprintf(fid,'%s%s',Eco(1,s).str{k,1},' ');
end
fprintf(fid,'\n');
else
fprintf(fid,'%s\n','Function Unknown');
end
end
end