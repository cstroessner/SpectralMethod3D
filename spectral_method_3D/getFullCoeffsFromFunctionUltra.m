function coeffs = getFullCoeffsFromFunctionUltra(f,n,lambda) 
% returns 3D Ultraspherical Coeffs C^{lambda} for f with resolution n
cheb = @(i,n) cos((i-1).*pi/(n-1));
ff = @(i,j,k) f(cheb(i,n(1)),cheb(j,n(2)),cheb(k,n(3)));
x = zeros([n(1),1,1]);
x(:,1,1) = 1:n(1);
X = repmat(x,1,n(2),n(3));
y = zeros([1,n(2),1]);
y(1,:,1) = 1:n(2);
Y = repmat(y,n(1),1,n(3));
z = zeros([1,1,n(3)]);
z(1,1,:) = 1:n(3);
Z = repmat(z,n(1),n(2),1);
coeffs = ff(X,Y,Z);
S1 = eye(n(1)); S2 = eye(n(2)); S3 = eye(n(3));
for iter = 1:lambda(1)
    S1 = getSUltra(n(1),iter-1) * S1;
end
for iter = 1:lambda(2)
    S2 = getSUltra(n(2),iter-1) * S2;
end
for iter = 1:lambda(3)
    S3 = getSUltra(n(3),iter-1) * S3;
end
coeffs = reshape(S1*Vals2ChebCoeffsMat(n(1))*reshape(coeffs,[n(1),n(2)*n(3)]),[n(1),n(2),n(3)]);
coeffs = permute(reshape(S2*Vals2ChebCoeffsMat(n(2))*reshape(permute(coeffs,[2,1,3]),[n(2),n(1)*n(3)]),[n(2),n(1),n(3)]),[2,1,3]);
coeffs = permute(reshape(S3*Vals2ChebCoeffsMat(n(3))*reshape(permute(coeffs,[3,2,1]),[n(3),n(2)*n(1)]),[n(3),n(2),n(1)]),[3,2,1]);
end