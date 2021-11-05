function F = Vals2ChebCoeffsMat(n)
% maps function evaluations at n Chebyshev nodes to Chebyshev coefficients
if n < 2
    warning('n too small')
end
F = zeros(n);
cheb = @(i,n) cos((i-1).*pi/(n-1));
xx = cheb(1:n,n);
T = @(i,x) cos((i-1)*acos(x));
for i = 1:n
    
    F(i,:) = T(i,xx);
    
end
F(:,1) = F(:,1)/2;
F(1,:) = F(1,:)/2;
F(:,n) = F(:,n)/2;
F(n,:) = F(n,:)/2;
F = (2/(n-1)).*F;
end