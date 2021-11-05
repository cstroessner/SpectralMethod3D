function D = getDUltra(n,lambda)
% Differentiation matrix mapping Chebyshev Coefficient to the C^{lambda}
% coefficients of the lamda-th derivative. See D_{lambda} in Townsend and
% Olver 2015.
if lambda == 0
    D = speye(n);
    return
end
D = spdiags((1:n)'-1,lambda,n,n); 
D = D.*(2^(lambda-1).*factorial(lambda-1));
end
