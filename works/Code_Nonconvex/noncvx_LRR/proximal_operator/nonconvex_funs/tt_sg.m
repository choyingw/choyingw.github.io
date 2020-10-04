function y = tt_sg(x, a,lambda)
% supergradient of lp penalty

x = abs(x) ;
epsilon = 0 ;

y = lambda*exp(-x.^3/a^3)/(a*sqrt(pi))/100+2;
    
%y = lambda*a*(x+epsilon).^(a-1) ; % 
% y = lambda*(x+epsilon).^(p-1) ;
