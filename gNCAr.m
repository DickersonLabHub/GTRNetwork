function [TFA CS]= gNCAr(Data)

E=Data.R;
A0=Data.C;

% Variable initialization
[N,M]= size(E);
L= size(A0,2);

if (size(A0,1) ~= N)
    error('The matrices E and A0 must have the same number of rows.');
end

A=A0;
A=A.*rand(size(A0));

P0=rand(L,M);
P=P0;

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

if (nz < L*(L-1))
    fprintf(1,'\n\nWarning: # of zeros in alpha < L*(L-1)\n');
    error('The system specified by ''A0'' and ''P0'' is NOT identifiable. Abort.');
elseif (R < L*(L-1))
        fprintf(1,'\n\nWarning: rank(B) < # L*(L-1)\n');
        error('The system specified by ''A0'' and ''P0'' is NOT identifiable. Abort.');
end

% Parameters and constraints for the EM algorithm
epsl= 1.0e-5;           % --> convergence threshold
stopeval=200;         %stop evaluation afer iteration >stopeval
% Identify the sparsity pattern of 'alpha'
sp_patternA = sign(A0);
sp_patternP = sign(P0);

% EM algorithm - main loop
delta_ssr_rel= 1.0;
ssr= 1e6; 
ssr_old =0;
iter = 1;
ss=[];

%Remove the zero columns before whitening data
zero_col=any(E);
zero_ind=find(zero_col);
E=E(:,zero_ind); 
M1=size(E,2);
P=P(:,zero_ind);


while ((abs(ssr_old -ssr) > epsl) && (iter < stopeval))
    % Step 1 - Re-estimation of 'P'
    for k = 1:M1
        index_P = find(sp_patternP(:,k) ~=0);
        A_red = A(:,index_P);
        P(index_P,k)= sol(A_red,E(:,k));
    end

    % Step 2 - Re-estimation of 'A' solving a sequence of least squares problems
    for n=1:N,
        % Build the low-dimension least square problem for the n_th gene
        index_A = find(sp_patternA(n,:) ~=0);
        P_red= P(index_A,:)';
        cP_red=cond(P_red);
        if cP_red > 100
            fprintf('\nGene %g is regulated by dependent TFs',n);
        end
        b=E(n,:)';
        alpha_red = P_red\b;
        A(n,index_A)= alpha_red';

    end
    iter = iter + 1;
    % (c) - Evaluate the current mismatch
    err = E-A*P;
    ssr_old = ssr;
    ssr = sum(abs(err(:)));
    ss=[ss;ssr];
    delta_ssr_rel= abs(ssr_old - ssr)/ssr_old;
    if 0==mod(iter,1000), 
        fprintf(1,'ssr= %.8f, delta_ssr_rel=%.10f\n',ssr,delta_ssr_rel);
    end
end

%Check the third criteria for solution P 
r=rank(P);
prob_g=[];
if r<L
    nullP=null(P');
    for n=1:N
        index_A = find(sp_patternA(n,:) == 0);   %row vector
        temp=nullP(index_A',:);
        temp_r=rank(temp);
        if temp_r<L-r
            prob_g=[prob_g;n];
        else
            cond_=cond(temp);
            if cond(temp)>10^7
                prob_g=[prob_g;n];
            end
        end
    end
end

%Add them back at the end
nzero_col=M-size(E,2);  %M is the original column number
if nzero_col
    Ptemp=P;
    P=[];
    count_=0;
    for j=1:M
        if zero_col(j)
            P(:,j)=Ptemp(:,count_+1);
            count_=count_+1;
        else
            P(:,j)=zeros(L,1);
        end
    end
end

if isempty(prob_g)
    fprintf(1,'\nThe solution satisfies the 3rd criterion\n');
else
    fprintf(1,'\nThe solution does not satisfy the 3rd criterion\n');
end
TFA=P;
CS=A;

%--------------------------------------------------------------------

function x = sol(A,b)
%solve least square problem by Tikhonov regularization method

[U,s,V]=csvd(A);  
lambda = l_curve(U,s,b);   
x = tikhonov(U,s,V,b,lambda); 


