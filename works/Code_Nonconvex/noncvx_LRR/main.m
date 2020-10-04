addpath(genpath(cd))
% Using functions from Lib-ADMM

% Optimization Parameter
opts.tol = 1e-9; 
opts.max_iter = 100;
opts.rho = 1.2;
opts.mu = 1e-3;
opts.max_mu = 1e9;
opts.DEBUG = 1;
opts.loss = 'l21';

%% LRR

lambda = 0.1;
T=2; % Number of samples
Ac1=zeros(1,T); Ac2=zeros(1,T); Ac3=zeros(1,T); Ac4=zeros(1,T); t1=zeros(1,1); t2=zeros(1,1); t3=zeros(1,1); t4=zeros(1,1);
par = [12 40 60]; % Parameter(p1,p2,p3) for proposed nonconvex method 
for t=1:T

% Generating the data of independent subspaces
d=10; amb=100; num_sub=10;S=zeros(100,100); idx=zeros(1,120);
for i=1:num_sub
    Sub=zeros(100,1); 
    r1=randperm(100,10); r2=rand(1,10);
    Sub(r1)=r2;
    Suba=[];
    for j=1:100
       Suba(:,end+1)=sprand(Sub);
    end
    SubS=Suba;
    [r,c]=find(SubS); C=unique(c); C=C(1:10);
    S(:,(i-1)*d+1:i*d)=SubS(:,C);
    idx(1,(i-1)*d+1:i*d)=i;
end


tic
[X1,E1,obj,err,iter,Err] = lrr(S,S,lambda,opts); % LRR
t1=t1+toc;
tic
[X2,E2,obj2,err2,iter2,Err2] = lrr2(S,S,lambda,opts); % LRR + dual moemntum
t2=t2+toc;
tic
[X3,E3,obj3,err3,iter3,Err3] = lrr3(S,S,lambda,opts,par); % nonconvex LRR
t3=t3+toc;
tic
[X4,E4,obj4,err4,iter4,Err4] = lrr4(S,S,lambda,opts,par); % nonconvex LRR + dual momentum
t4=t4+toc;

% Calculating the clusering accuracy

[W1,Wo1,acc1]=showW2(X1,idx); Ac1(t)=acc1; 
[W2,Wo2,acc2]=showW2(X2,idx); Ac2(t)=acc2;
[W3,Wo3,acc3]=showW2(X3,idx); Ac3(t)=acc3;
[W4,Wo4,acc4]=showW2(X4,idx); Ac4(t)=acc4;
end

As1=mean(Ac1)
t1=t1/T
As2=mean(Ac2)
t2=t2/T
As3=mean(Ac3)
t3=t3/T
As4=mean(Ac4) 
t4=t4/T



