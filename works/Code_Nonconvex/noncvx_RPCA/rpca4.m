function [L,S,obj,err,iter,Err] = rpca4(X,lambda,opts,ori)

% nonconvexRPCA + dual momentum adapted from LibADMM

% Solve the Robust Principal Component Analysis minimization problem by M-ADMM
%
% min_{L,S} ||L||_*+lambda*loss(S), s.t. X=L+S
% loss(S) = ||S||_1 or ||S||_{2,1}
%
% ---------------------------------------------
% Input:
%       X       -    d*n matrix
%       lambda  -    >0, parameter
%       opts    -    Structure value in Matlab. The fields are
%           opts.loss       -   'l1' (default): loss(S) = ||S||_1 
%                               'l21': loss(S) = ||S||_{2,1}
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.DEBUG      -   0 or 1
%
% Output:
%       L       -    d*n matrix
%       S       -    d*n matrix
%       obj     -    objective function value
%       err     -    residual 
%       iter    -    number of iterations


tol = 1e-8; 
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;
DEBUG = 0;
loss = 'l1';

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'loss');        loss = opts.loss;            end
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end


[d,n] = size(X);

L = zeros(d,n);
S = L;
Y = L;
Err=[];
iter = 0; alpha=1; S_hat = S; Y_hat = Y; c = 10; 
for iter = 1 : max_iter
    Lk = L;
    Sk = S;
    % update L
    [L,nuclearnormL] = prox_nuclear2(-S+X-Y_hat/mu,1/mu,[10 40 70]);
    % update S
    if strcmp(loss,'l1')
        S = prox_l1(-L+X-Y_hat/mu,lambda/mu);
    elseif strcmp(loss,'l21')
        S = prox_l21(-L+X-Y_hat/mu,lambda/mu);
    else
        error('not supported loss function');
    end
  
    dY = L+S-X;
    chgL = max(max(abs(Lk-L)));
    chgS = max(max(abs(Sk-S)));
    chg = max([chgL chgS max(abs(dY(:)))]);
    err = norm(dY,'fro');
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            obj = nuclearnormL+lambda*comp_loss(S,loss);
            err = norm(dY,'fro');
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err) ', c =' num2str(c)]); 
        end
    end
    
    if err < tol
        break;
    end 
    Yk = Y;
    Y = Y_hat + mu*dY;
    
    alphak = alpha;
    alpha = (1+sqrt(1+4*alpha^2))/2;
    Y_hat = Y+((alphak-1)/alpha)*(Y-Yk);
    mu = min(rho*mu,max_mu);
    Err(end+1)= norm(ori-L,'fro');
end
obj = nuclearnormL+lambda*comp_loss(S,loss);
err = norm(dY,'fro');

function out = comp_loss(E,loss)

switch loss
    case 'l1'
        out = norm(E(:),1);
    case 'l21'
        out = 0;
        for i = 1 : size(E,2)
            out = out + norm(E(:,i));
        end
end