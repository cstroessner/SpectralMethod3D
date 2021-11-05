% clear
% n = [10,10,10]
% 
% uCoeffs = rand(n);
% uTrue = @(x,y,z)  0.*x.*y.*z;
% rightbc = @(y,z) uTrue(1,y,z);
% leftbc = @(y,z) uTrue(-1,y,z);
% topbc = @(x,z) uTrue(x,1,z);
% bottombc = @(x,z) uTrue(x,-1,z);
% frontbc = @(x,y) uTrue(x,y,-1);
% backbc = @(x,y) uTrue(x,y,1);
% 
% [T1,F1,T2,F2,T3,F3] = getBoundaryConditionMatrices(n,rightbc,leftbc,topbc,bottombc,frontbc,backbc);
% R1 = max(abs(tprod(uCoeffs,T1,eye(n(2)),eye(n(3)))-F1),[],'all');
% R2 = max(abs(tprod(uCoeffs,eye(n(1)),T2,eye(n(3)))-F2),[],'all');
% R3 = max(abs(tprod(uCoeffs,eye(n(1)),eye(n(2)),T3)-F3),[],'all');
% ResidualNorm = max([R1,R2,R3])
% 
% uCoeffs = projectOntoDirichlet1(uCoeffs,rightbc,leftbc,topbc,bottombc,frontbc,backbc);
% 
% R1 = max(abs(tprod(uCoeffs,T1,eye(n(2)),eye(n(3)))-F1),[],'all');
% R2 = max(abs(tprod(uCoeffs,eye(n(1)),T2,eye(n(3)))-F2),[],'all');
% R3 = max(abs(tprod(uCoeffs,eye(n(1)),eye(n(2)),T3)-F3),[],'all');
% ResidualNorm = max([R1,R2,R3])

function uCoeffs = projectOntoDirichlet(uCoeffsUltra,rightbc,leftbc,topbc,bottombc,frontbc,backbc)
% TODO use FFT instead of F matrix

% uuCoeffs = uCoeffs;

n = size(uCoeffsUltra);

% get Cheby instead of Ultra coeffs
uCoeffs = invtprod(uCoeffsUltra,getSUltra(n(1),1),getSUltra(n(2),1),getSUltra(n(3),1));
uCoeffs = invtprod(uCoeffs,getSUltra(n(1),0),getSUltra(n(2),0),getSUltra(n(3),0));


% get fun values at cheby nodes
F1 = Vals2ChebCoeffsMat(n(1));
F2 = Vals2ChebCoeffsMat(n(2));
F3 = Vals2ChebCoeffsMat(n(3));
uVals = invtprod(uCoeffs,F1,F2,F3);

% uuVals = uVals;

% fix boundaries 
cheb = @(i,n) cos((i-1).*pi/(n-1));

for k = 1:n(3)
uVals(1,:,k) = leftbc(cheb(1:n(2),n(2)),cheb(k,n(3)));
uVals(end,:,k) = rightbc(cheb(1:n(2),n(2)),cheb(k,n(3)));
uVals(:,1,k) = bottombc(cheb(1:n(1),n(1)),cheb(k,n(3)));
uVals(:,end,k) = topbc(cheb(1:n(1),n(1)),cheb(k,n(3)));
end
for j = 1:n(1)
    uVals(:,j,1) = frontbc(cheb(1:n(1),n(1)),cheb(j,n(2)));
    uVals(:,j,end) = backbc(cheb(1:n(1),n(1)),cheb(j,n(2)));
end

% project back
uCoeffs = tprod(uVals,F1,F2,F3);

% diffVals = max(abs(uVals-uuVals),[],'all')
% diffCoeffs = max(abs(uCoeffs-uuCoeffs),[],'all')
end

