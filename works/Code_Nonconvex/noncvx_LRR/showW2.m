function [L,L2,acc]=showW2(X,label)

% Adapted from LibADMM and use clustering accuracy from  Praisan Padungweang

gnd = label;
K = max(gnd);

%post processing
[U,S,V] = svd(X,'econ');
S = diag(S);
r = sum(S>1e-4*S(1));
U = U(:,1:r);S = S(1:r);
U = U*diag(sqrt(S));
U = normr(U);
L = (U*U').^4;
L2=L;

D = diag(1./sqrt(sum(L,2)));
L = D*L*D;
[U,S,V] = svd(L);
V = U(:,1:K);
V = D*V;
idx = kmeans(V,K,'emptyaction','singleton','replicates',100,'display','off');
acc = AccMeasure(idx,gnd)/100;
disp(['seg acc=' num2str(acc)]);