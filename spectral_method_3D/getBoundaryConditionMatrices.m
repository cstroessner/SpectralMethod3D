function [T1,F1,T2,F2,T3,F3] = getBoundaryConditionMatrices(n,rightbc,leftbc,topbc,bottombc,frontbc,backbc)
% output boundary conditions in the form X x_1 T_1 = F_1, X x_2 T_2 = F_2,
% X x_3 T_3 = F_3 for the given boundary functions
% the matrices T have identity as principal submatrices

% 1 (left + right)
coeffsRight = getFull2DCoeffsFromFunction(rightbc,[n(2),n(3)]);
coeffsLeft = getFull2DCoeffsFromFunction(leftbc,[n(2),n(3)]);
F1 = zeros([2,n(2),n(3)]);
F1(1,:,:) = coeffsRight;
F1(2,:,:) = coeffsLeft;
T1 = [evalChebyPolys(n(1),1)';evalChebyPolys(n(1),-1)'];
% 2 (top + bottom)
coeffsTop = getFull2DCoeffsFromFunction(topbc,[n(1),n(3)]);
coeffsBot = getFull2DCoeffsFromFunction(bottombc,[n(1),n(3)]);
F2 = zeros([n(1),2,n(3)]);
F2(:,1,:) = coeffsTop;
F2(:,2,:) = coeffsBot;
T2 = [evalChebyPolys(n(2),1)';evalChebyPolys(n(2),-1)'];
% 3 (front + back)
coeffsFront = getFull2DCoeffsFromFunction(frontbc,[n(1),n(2)]);
coeffsBack = getFull2DCoeffsFromFunction(backbc,[n(1),n(2)]);
F3 = zeros([n(1),n(2),2]);
F3(:,:,1) = coeffsBack;
F3(:,:,2) = coeffsFront;
T3 = [evalChebyPolys(n(3),1)';evalChebyPolys(n(3),-1)'];
% modify such that principal submatrices in T are identities.
T1princsub = T1(1:2,1:2);
T1 = T1princsub\T1;
F1 = reshape( T1princsub\reshape(F1,[2,n(2)*n(3)]),[2,n(2),n(3)]);
T2princsub = T2(1:2,1:2);
T2 = T2princsub\T2;
F2 = permute(reshape(T1princsub\reshape(permute(F2,[2,1,3]),[2,n(1)*n(3)]),[2,n(1),n(3)]),[2,1,3]);
T3princsub = T3(1:2,1:2);
T3 = T3princsub\T3;
F3 = permute(reshape(T3princsub\reshape(permute(F3,[3,2,1]),[2,n(2)*n(1)]),[2,n(2),n(1)]),[3,2,1]);
end