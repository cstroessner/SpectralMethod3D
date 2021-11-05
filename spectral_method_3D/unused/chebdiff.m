function [D,x] = chebdiff(n)
% maps function evaluations at Chebyshev nodes to function evaluations of
% the derivative at the same Chebyshev nodes
% taken from https://people.maths.ox.ac.uk/trefethen/spectral.html 
% see Lloyd N. Trefethen, Spectral Methods in MATLAB, SIAM, Philadelphia, 2000
n = n-1;
if n==0, D=0; x=1; return, end
x = cos(pi*(0:n)/n)';
c = [2; ones(n-1,1); 2].*(-1).^(0:n)';
X = repmat(x,1,n+1);
dX = X-X';
D  = (c*(1./c)')./(dX+(eye(n+1)));
D  = D - diag(sum(D,2));
end