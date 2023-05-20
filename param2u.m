function u = param2u( qNo, alpha, tau)
%
% u = param2u( F0, Fm, p, sigma, PAR, alpha, tau )
% transform parameters to the fit variable u
%
 
 
nExp = length(tau);
 
nw = 1; 
u = zeros( nw + 2 * nExp -1, 1 ); 
u(1) = qNo; 

 
tau = tau(:);
[tau,ind] = sort(tau,'ascend'); 
alpha = alpha(ind); 
 
u(nw+(1:nExp)) = [tau(1); diff(tau) ];  
 
if nExp > 1
    alpha = alpha / sum(alpha);  % enforce sum(alpha)=1
    i0 = nw+nExp;
    u(i0+1) = alpha(1);
    for n=2:nExp-1
        % u(i0+n) is alpha(n) as a fraction of alpha budget left
        % after alpha(1),...,alpha(n-1)
        budget = 1 - sum(alpha(1:n-1));
        if budget > 0
            u(i0+n) = alpha(n) /  budget;
        else
            u(i0+n) = 0;
        end
    end
end
 
end


