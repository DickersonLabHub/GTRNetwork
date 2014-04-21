function gNCA_uniq = checkgNCA( Expression,Regulatory )
%   Check gNCA solution 
%   by Fu, Yao(Al)
Data0=MatchGenes(Expression,Regulatory);
Data=Subnetwork(Data0);

E=Data.R;
A0=Data.C;
[N,M]= size(E);
L= size(A0,2);

if (size(A0,1) ~= N)
    error('The matrices E and A0 must have the same number of rows.');
end

P0=rand(L,M);
%Check the rank condition
R= 0;   % --> initialize rank of the augmented sparse system
nz= 0;  % --> initialize total # of zeros in the augmented system

if rank(A0)<L, error('A0 is not full column rank'),
end

for l = 1:L
    B=[];
    %Construct reduced matrix P
    pr = [];
    for i = 1:L
        nr = 0;
        if i ~= l
            r = find(P0(i,:) == 0);
            nr = length(r);
            pl = zeros(nr,L);
            pl(:,i)=P0(l,r)';
            nz=nz+nr;
            pr=[pr;pl];    
        end;         
    end;

    v= find(A0(:,l)==0);
    B=[A0(v,:);pr];
    nz= nz+length(v);
    R= R+ rank(B);
end

gNCA_uniq=1;
if (nz < L*(L-1))
    fprintf(1,'\n\nWarning: # of zeros in alpha < L*(L-1)\n');
    gNCA_uniq=0;
elseif (R < L*(L-1))
        fprintf(1,'\n\nWarning: rank(B) < # L*(L-1)\n');
        gNCA_uniq=0;
end


end

