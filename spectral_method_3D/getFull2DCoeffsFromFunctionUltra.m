function coeffs = getFull2DCoeffsFromFunctionUltra(f,n,lambda)
% computed 2D Chebyshev coeffs for f with resolution n (used for evaluating the 2Dboudary functions)
cheb = @(i,n) cos((i-1).*pi/(n-1));
ff = @(i,j,k) f(cheb(i,n(1)),cheb(j,n(2)));
x = zeros([n(1),1,1]);
x(:,1,1) = 1:n(1);
X = repmat(x,1,n(2));
y = zeros([1,n(2),1]);
y(1,:,1) = 1:n(2);
Y = repmat(y,n(1),1);
coeffs = ff(X,Y);

S1 = eye(n(1)); S2 = eye(n(2));
for iter = 1:lambda(1)
    S1 = getSUltra(n(1),iter-1) * S1;
end
for iter = 1:lambda(2)
    S2 = getSUltra(n(2),iter-1) * S2;
end

coeffs = S1*Vals2ChebCoeffsMat(n(1))*coeffs;
coeffs = permute(S2*Vals2ChebCoeffsMat(n(2))*reshape(permute(coeffs,[2,1]),[n(2),n(1)]),[2,1]);
end