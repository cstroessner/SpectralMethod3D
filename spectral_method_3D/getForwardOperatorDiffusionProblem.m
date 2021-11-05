function [opLCP,lambda] = getForwardOperatorDiffusionProblem(a1,a2,a3,n)
% Computes the forward differential operator as map from 3D Chebyshev
% coefficients to 3D C^{lambda} ultraspherical coefficients. The returned matrices need to
% be multiplied in the corresponding modes and summed up afterwards. See
% applyForwardOperator for details.

lambda = [2,2,2]; 

% get cheby coeffs of 1D functions
cheb = @(i,n) cos((i-1).*pi/(n-1));
aa1 = @(i) a1(cheb(i,n(1)));
x = 1:n(1);
a1Coeffs = Vals2ChebCoeffsMat(n(1))*aa1(x)';

cheb = @(i,n) cos((i-1).*pi/(n-1));
aa2 = @(i) a2(cheb(i,n(2)));
x = 1:n(2);
a2Coeffs = Vals2ChebCoeffsMat(n(2))*aa2(x)';

cheb = @(i,n) cos((i-1).*pi/(n-1));
aa3 = @(i) a3(cheb(i,n(3)));
x = 1:n(3);
a3Coeffs = Vals2ChebCoeffsMat(n(3))*aa3(x)';

% get matrices for the operator
opLCP{1,1} = -getSUltra(n(1),1)*getDUltra(n(1),1)*(getSUltra(n(1),0)\(MultiplicationMatrix(a1Coeffs,1)*getDUltra(n(1),1))); %TODO it would be more elegant to use a differentiation matrix form C1 to C2
opLCP{2,1} = MultiplicationMatrix(a2Coeffs,2)*getSUltra(n(2),1)*getSUltra(n(2),0);
opLCP{3,1} = MultiplicationMatrix(a3Coeffs,2)*getSUltra(n(3),1)*getSUltra(n(3),0);

opLCP{1,2} = MultiplicationMatrix(a1Coeffs,2)*getSUltra(n(1),1)*getSUltra(n(1),0);
opLCP{2,2} = -getSUltra(n(2),1)*getDUltra(n(2),1)*(getSUltra(n(2),0)\(MultiplicationMatrix(a2Coeffs,1)*getDUltra(n(2),1)));
opLCP{3,2} = MultiplicationMatrix(a3Coeffs,2)*getSUltra(n(3),1)*getSUltra(n(3),0);

opLCP{1,3} = MultiplicationMatrix(a1Coeffs,2)*getSUltra(n(1),1)*getSUltra(n(1),0);
opLCP{2,3} = MultiplicationMatrix(a2Coeffs,2)*getSUltra(n(2),1)*getSUltra(n(2),0);
opLCP{3,3} = -getSUltra(n(3),1)*getDUltra(n(3),1)*(getSUltra(n(3),0)\(MultiplicationMatrix(a3Coeffs,1)*getDUltra(n(3),1)));
end