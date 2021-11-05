function [T1,F1,T2,F2,T3,F3] = getBoundaryConditionMatricesOriginal(n,rightbc,leftbc,topbc,bottombc,frontbc,backbc)
% output boundary conditions in the form X x_1 T_1 = F_1, X x_2 T_2 = F_2,
% X x_3 T_3 = F_3 for the given boundary functions
% the matrices T do not have identity as principal submatrices

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
end