% test on the matrix completion problem
% Adapted from  C. Lu, J. Tang, S. Yan, and Z. Lin, ¡§Nonconvex nonsmooth low rank minimization via iteratively reweighted nuclear norm,¡¨ IEEE Trans. Image Process., vol. 25, no. 2, pp. 829-839, Feb. 2016.
currpath = cd ;
addpath(genpath(currpath)) ; 
testNum=100; % Number of testing data
err=zeros(10,13,testNum); times=err; % Performance evaluation idx:(rank,method,testNum);


parfor T=1:testNum % index of testing number 
newDataFlag = 1 ;
rrr=10:2:28; % rank of underlying matrix

for kk=1:10 % index of rrr

if newDataFlag == 1    

    m = 150 ; % matrix size
    n = m ;
    r = rrr(kk);
    ML = (randn(m,r)); MR = (randn(n,r));
    p = 0.5 * m* n ;
    rho_s = p / (m * n);    
    [I J col omega] = myRandsample(m, n, p);
    V = UVtOmega(ML, MR, I, J, col);    
    D = spconvert([I,J,V; m,n,0]);
%     clear I J col;    
end

IDX = omega ;
sizeX = [m,n] ;
M = opRestriction(prod(sizeX), IDX);
X = ML*MR' ; % Low-Rank matrix
x = X(:);
y = M(x,1); % Corruption

% Parameter for all the methods

tstD=reshape(M(y,2),[m,n]);
[trow, tcol, tval] = find(tstD);

para=[];
para.test.row = trow;
para.test.col = tcol;
para.test.data = tval;
para.test.m = size(tstD, 1);
para.test.n = size(tstD, 2);
para.tol = 1e-8;
para.maxIter = 1000;
para.tau = 1.01;
para.decay = 0.9;
lambda=0.1;


for i=6:6 % Choosing the methods i=1~5 are other nonconvex penalization i=6 is proposed nonconvex surrogate i=7~13 are ohter methods
if i==1
 fun = 'lp' ;        gamma = 0.5 ;
elseif i==2
 fun = 'scad' ;      gamma = 100 ;
elseif i==3
 fun = 'mcp' ; gamma = 10 ;
elseif i==4
 fun = 'cappedl1' ; gamma = 100 ;
elseif i==5
 fun = 'etp' ;  gamma = 0.1 ;
elseif i==6
fun = 'linear'; gamma=[5,50,60];
end

if T~=8 && T~=9
if i==7 || i==8
EOT=1;  
if i==8
    EOT=2;
end
xx=X; I_miss=reshape(M(M(xx(:),1),2),[m,n]); K=false(m*n,1); K(omega)=true; Omega=reshape(K,[m,n]);
tol=  10^-6; U_BL =   ones(1,1)*10^-3;  ruo =1.05 ; L=20; X_end=1.5; tr_jd=12; Itmax=100;
tic;
[ETNNR_M1,num1,ETNNR_M2]=ETNNR(xx,I_miss,Omega,  Itmax,tol,ruo, U_BL, L,X_end, tr_jd, EOT,4);
times(kk,i,T) = times(kk,i,T)+toc;
XRec= ETNNR_M1(:,:,:,end);
%XRec = max(0,XRec); XRec = min(1,XRec);
end
end

if i==9
tic;
[ U, S, V, ~ ] = mc_alt( tstD, lambda, para );
times(kk,i,T) = times(kk,i,T)+toc;
XRec = U*S*V';
end

if i==10
tic;
[U, S, V,~] = APGMatComp(tstD,lambda, para );
times(kk,i,T) = times(kk,i,T)+toc;
XRec = U*S*V';
end

if i==11
para.maxR = 80;
tic;
para.speedup = 1;
para.exact = 1;
[ U, S, V, ~ ] = SoftImpute(tstD,lambda, para);
times(kk,i,T) = times(kk,i,T)+toc;
XRec = U*S*V';
end

if i==12
K=false(m*n,1); K(omega)=true; Omega=reshape(K,[m,n]);   lower_R = 9; upper_R = 9;  
tic;
[admmret, ~]= admm_pic('test',X,Omega,0,lower_R,upper_R);
times(kk,i,T) = times(kk,i,T)+toc;
XRec = admmret.recover;
end

if i==13
K=false(m*n,1); K(omega)=true; Omega=reshape(K,[m,n]);   lower_R = 9; upper_R = 9;  
tic;
[apglret, ~]= apgl_pic('test',X,Omega,0,lower_R,upper_R);
times(kk,i,T) = times(kk,i,T)+toc;
XRec = apglret.recover;
end

if i<7
tic;
[XRec,~,~] = IRNN(fun,y,M,m,n,gamma,X) ;
times(kk,i,T) = times(kk,i,T)+toc;
end

err(kk,i,T)=  norm(XRec-X,'fro')/norm(X,'fro');

end
end
end

averageErr=mean(err,3);
averageTime=mean(times,3);