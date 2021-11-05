function coeffs = getEvalsFromFunction(f,n)
% returns 3D Chebyshev point evals for f with resolution n
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
end