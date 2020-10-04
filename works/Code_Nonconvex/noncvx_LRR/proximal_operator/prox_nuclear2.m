function [X,nuclearnorm] = prox_nuclear2(B,lambda,par)

% Proximal operator from Lib-ADMM with adaption using nonconvex penalization

% The proximal operator of the nuclear norm of a matrix
% 
% min_X lambda*||X||_*+0.5*||X-B||_F^2


[U,W,V] = MySVD(B);
W = diag(W);
w = linear_sg(W,par,lambda); %12 5 10
%w = etp_sg(W,0.01,lambda);
%w = mcp_sg(W,10,lambda);
W = W - w;
svp=length(find(W>0));
X = U(:,1:svp)*diag(W(1:svp))*V(:,1:svp)';
nuclearnorm = sum(W);
if svp<1
    X = zeros(size(B));
    nuclearnorm = 0;
end
