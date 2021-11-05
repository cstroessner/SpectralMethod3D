function coeffs = get1DCoeffsFromFunctionUltra(f,n,lambda)
% computes C-lambda ultraspherical coefficients for one-dimensional
% functions
cheb = @(i,n) cos((i-1).*pi/(n-1));
evals = @(i) f(cheb(i,n));
coeffs = Vals2ChebCoeffsMat(n)*evals(1:n)';
for k = 0:lambda-1
   coeffs = getSUltra(n,k)*coeffs; 
end
end

