% f = @(x,y,z) (x-1).*(x+1).*(y-1).*(y+1).*(z-1).*(z+1)+1;
% 
% bcFun = @(x,y,z) 3.*x.*y.*z+sin(z);
% rightbc = @(y,z) bcFun(1,y,z);
% leftbc = @(y,z) bcFun(-1,y,z);
% topbc = @(x,z) bcFun(x,1,z);
% bottombc = @(x,z) bcFun(x,-1,z);
% frontbc = @(x,y) bcFun(x,y,-1);
% backbc = @(x,y) bcFun(x,y,1);
% 
% 
% uCoeffsUlt = getFullCoeffsFromFunctionUltra(f,[10,10,10],[2,2,2]);
% uCoeffs = projectDirichletBCFunctions1(uCoeffsUlt,rightbc,leftbc,topbc,bottombc,frontbc,backbc);
% 
% sanityCheck = max(abs([funeval(uCoeffs,0,0,1)-bcFun(0,0,1),funeval(uCoeffs,0,-1,0)-bcFun(0,-1,0),funeval(uCoeffs,1,0,1)-bcFun(1,0,1)]))
% 

function uCoeffs = projectDirichletBCFunctions(uCoeffs,rightbc,leftbc,topbc,bottombc,frontbc,backbc)
% projects a function given in terms of C 2 ultraspherical coefficients onto a
% function with given boundary conditions
% afterwards the coefficients are mapped to chebyshev coefficients

% get evaluation matrices
lambda = 2; k = 2;
n = size(uCoeffs);
C1 = evalUltraPolynomials(n(1),lambda,[1,-1]);
Cprincsub1 = C1(end-1:end,end-1:end);
C1 = Cprincsub1\C1;
C2 = evalUltraPolynomials(n(2),lambda,[1,-1]);
Cprincsub2 = C2(end-1:end,end-1:end);
C2 = Cprincsub2\C2;
C3 = evalUltraPolynomials(n(3),lambda,[1,-1]);
Cprincsub3 = C3(end-1:end,end-1:end);
C3 = Cprincsub3\C3;
T1 = C1;
T2 = C2;
T3 = C3;

% get desired boundary values
coeffsRight = getFull2DCoeffsFromFunctionUltra(rightbc,[n(2),n(3)],[2,2]);
coeffsLeft = getFull2DCoeffsFromFunctionUltra(leftbc,[n(2),n(3)],[2,2]);
F1 = zeros([2,n(2),n(3)]);
F1(1,:,:) = coeffsRight;
F1(2,:,:) = coeffsLeft;
coeffsTop = getFull2DCoeffsFromFunctionUltra(topbc,[n(1),n(3)],[2,2]);
coeffsBot = getFull2DCoeffsFromFunctionUltra(bottombc,[n(1),n(3)],[2,2]);
F2 = zeros([n(1),2,n(3)]);
F2(:,1,:) = coeffsTop;
F2(:,2,:) = coeffsBot;
coeffsFront = getFull2DCoeffsFromFunctionUltra(frontbc,[n(1),n(2)],[2,2]);
coeffsBack = getFull2DCoeffsFromFunctionUltra(backbc,[n(1),n(2)],[2,2]);
F3 = zeros([n(1),n(2),2]);
F3(:,:,1) = coeffsBack;
F3(:,:,2) = coeffsFront;
F1 = reshape( Cprincsub1\reshape(F1,[2,n(2)*n(3)]),[2,n(2),n(3)]);
F2 = permute(reshape(Cprincsub2\reshape(permute(F2,[2,1,3]),[2,n(1)*n(3)]),[2,n(1),n(3)]),[2,1,3]);
F3 = permute(reshape(Cprincsub3\reshape(permute(F3,[3,2,1]),[2,n(2)*n(1)]),[2,n(2),n(1)]),[3,2,1]);

% perform the coefficient change
uCoeffs(end-1:end,1:end-2,1:end-2) = F1(:,1:end-2,1:end-2) - reshape(T1(:,1:end-2)*reshape(uCoeffs(1:end-2,1:end-2,1:end-2),[n(1)-k,(n(2)-k)*(n(3)-k)]),[2,n(2)-k,n(3)-k]);
uCoeffs(:,end-1:end,1:end-2) = F2(:,:,1:end-2) - permute(reshape(T2(:,1:end-2)*reshape(permute(uCoeffs(:,1:end-2,1:end-2),[2,1,3]),[n(2)-k,(n(1))*(n(3)-k)]),[2,n(1),n(3)-k]),[2,1,3]);
uCoeffs(:,:,end-1:end) = F3(:,:,:) - permute(reshape(T3(:,1:end-2)*reshape(permute(uCoeffs(:,:,1:end-2),[3,2,1]),[n(3)-k,(n(2))*(n(1))]),[2,n(2),n(1)]),[3,2,1]);

% C = evalUltraPolynomials(n(1),lambda,[1,-1]);
% err = max( [max(abs(tprod(uCoeffs,eye(10),C,eye(10))-F2),[],'all'), ...
%     max(abs(tprod(uCoeffs,eye(10),eye(10),C)-F3),[],'all'), ...
%     max(abs(tprod(uCoeffs,C,eye(10),eye(10))-F1),[],'all')])

% map to cheby coeffs 
uCoeffs = invtprod(uCoeffs,getSUltra(n(1),1),getSUltra(n(2),1),getSUltra(n(3),1));
uCoeffs = invtprod(uCoeffs,getSUltra(n(1),0),getSUltra(n(2),0),getSUltra(n(3),0));
end

