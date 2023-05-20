function  [qNo, alpha, tau] = u2param( u )
%
% [F0, Fm, p, sigma, PAR, alpha, tau] = u2param( u )
% transform the fit variable u to parameters 
%
 
 
qNo = u(1);
 
nw = 1; 
 
nExp = ( length(u)-nw +1) / 2; 
assert(nExp==round(nExp)); % something is wrong if nExp is not integer
 
% tau:  relaxation times
% This transformation makes tau an ascending sequence
% tau(3) >= tau(2) >= tau(1)
% Otherwise the optimization problem would be degenerate
% tau_diff(1) = tau(1); 
% tau_diff(2) = tau(2)-tau(1);
% tau_diff(3) = tau(3)-tau(2)
tau_diff = u(nw+(1:nExp)); 
tau = cumsum( tau_diff ); 
tau = reshape( tau, [1,nExp]); 
 
% alpha's are determined from v's
% 0 <= v <= 1
v = u(nw+1+nExp:end); 
assert(length(v)==nExp-1); 
 
% The following assures that: 
% 0 <= alpha <= 1 and sum(alpha)=1
%
% alpha(1) = v(1)
% 1-alpha(1) is teh budget left for the rest of the alphas
% v(2) is the fraction of that budget used for alpha(2)
% alpha(2) = v(2) * (1-alpha(1))
% This will work for any number of alphas
 
alpha = zeros( 1, nExp ); 
 
if length(v) >= 1
    alpha(1) = v(1);
end
for n=2:nExp-1
    alpha(n) = ( 1-sum(alpha(1:n-1)) ) * v(n); 
end
alpha(end) = 1-sum(alpha(1:end-1)); 
 
end
