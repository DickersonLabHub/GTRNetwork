function [ result ] = AssambleResult(result_folder,includingSigma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
d=dir(result_folder);
if includingSigma==0
    for i=1:4
        for k=1:5
            for j=1:4
            load([result_folder '/' d(20*(i-1)+2+5*(j-1)+k).name]);
            end
            Pnone=debug3(Pnone);
            Psvd=debug3(Psvd);
            Ppls=debug3(Ppls);
            Pem=debug3(Pem);
            result(i,k).None=SumResults(Pnone);
            result(i,k).SVD=SumResults(Psvd);
            result(i,k).PLS=SumResults(Ppls);
            result(i,k).EM=SumResults(Pem);
            result(i,k).note=[ 'Result for ininital network input setting at ' num2str(i*20+10) 'percent of total known regulatory links. Replication' num2str(k)];
        fprintf('%i%s\n',5*(k+(i-1)*5),'% complated');
        end
    end
else
    for i=1:4
        for k=1:5
            for j=1:3
            load([result_folder '/' d(15*(i-1)+2+5*(j-1)+k).name]);
            end
            Psvd=debug3(Psvd);
            Ppls=debug3(Ppls);
            Pem=debug3(Pem);
            result(i,k).SVD=SumResults(Psvd);
            result(i,k).PLS=SumResults(Ppls);
            result(i,k).EM=SumResults(Pem);
            result(i,k).note=[ 'Result for ininital network input setting at ' num2str(i*20+10) 'percent of total known regulatory links. Replication' num2str(k)];
        fprintf('%i%s\n',5*(k+(i-1)*5),'% complated');
        end
        
    end
end

