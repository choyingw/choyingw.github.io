addpath(genpath(cd))
% Using functions from Lib-ADMM
% Data initialization
d = 10;
na = 200;
nb = 100;

A = randn(d,na);
X = randn(na,nb);
B = A*X;
b = B(:,1);

%Optimization Parameter
opts.tol = 1e-9; 
opts.max_iter = 100;
opts.rho = 1.2;
opts.mu = 1e-3;
opts.max_mu = 1e9;
opts.DEBUG = 1;
opts.loss = 'l1'; 
opts.DEBUG = 1;

%% RPCA

M=20; % How many testing samples
Maxit=100; % Iterations
t1=0; t2=t1; t3=0; t4=0; Ek1=zeros(M,Maxit); Ek2=Ek1;Ek3=Ek1; Ek4=Ek1; 
n=100;
for k=1:M

n = 100; % Low-rank matrix generalization
r = 10;
X = rand(n,r)*rand(r,n);

Err = 0.01*rand(100,100);
K = rand(100,100); Err(K<=0.8)=0;
B = X+Err;
B2=X+sqrt(Err);
lambda = 1/sqrt(n);

tic
[L,S,obj,err,iter,Err] = rpca(B,lambda,opts,X); % RPCA
t1=t1+toc;
Err(end+1:Maxit)=0; 
Ek1(k,:)=Err;

tic
[L2,S2,obj2,err2,iter2,Err2] = rpca2(B,lambda,opts,X); %RPCA with dual momentum
t2=t2+toc;
Err2(end+1:Maxit)=0;
Ek2(k,:)=Err2;

tic
[L3,S3,obj3,err3,iter3,Err3] = rpca3(B,lambda,opts,X); % nonconvex RPCA
t3=t3+toc;
Err3(end+1:Maxit)=0;
Ek3(k,:)=Err3;

tic
[L4,S4,obj4,err4,iter4,Err4] = rpca4(B,lambda,opts,X); % nonconvex RPCA with dual momentum
t4=t4+toc;
Err4(end+1:Maxit)=0;
Ek4(k,:)=Err4;

end
t=[t1 t2 t3 t4]'/M;
Ek1=mean(Ek1,1); Ek2=mean(Ek2,1); Ek3=mean(Ek3,1); Ek4=mean(Ek4,1);
figure();
semilogy(1:Maxit,Ek1,1:Maxit,Ek2,1:Maxit,Ek3,1:Maxit,Ek4);
legend('RPCA','RPCA + dual momentum','nonconvex RPCA','nonconvex RPCA + dual momentum');

