function coeffs = get1DCoeffsFromFunction(f,n)
% computes Chebyshev coeffs for 1D functions
cheb = @(i,n) cos((i-1).*pi/(n-1));
evals = @(i) f(cheb(i,n));
coeffs = Vals2ChebCoeffsMat(n)*evals(1:n)';
end

