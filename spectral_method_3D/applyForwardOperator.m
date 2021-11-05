function LuCoeffs = applyForwardOperator(opLCP,uCoeffs)
% Applies the operators obtained from getForwardOperator to Chebyshev
% Coeffs. Depending on the matrices in opLCP we obtain 3D
% Chebyshev or Ultraspherical Coefficients.

n = [size(opLCP{1,1},1),size(opLCP{2,1},1),size(opLCP{3,1},1)];
r = size(opLCP,2);

LuCoeffs = zeros(size(uCoeffs));
for i = 1:r
    LiCoeffs = uCoeffs;
    LiCoeffs = reshape(opLCP{1,i}*reshape(LiCoeffs,[n(1),n(2)*n(3)]),[n(1),n(2),n(3)]);
    LiCoeffs = permute(reshape(opLCP{2,i}*reshape(permute(LiCoeffs,[2,1,3]),[n(2),n(1)*n(3)]),[n(2),n(1),n(3)]),[2,1,3]);
    LiCoeffs = permute(reshape(opLCP{3,i}*reshape(permute(LiCoeffs,[3,2,1]),[n(3),n(2)*n(1)]),[n(3),n(2),n(1)]),[3,2,1]);
    LuCoeffs = LuCoeffs + LiCoeffs;
end
end
