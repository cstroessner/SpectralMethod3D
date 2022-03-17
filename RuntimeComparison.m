% this generates the values for the table in Section 5.1

clear
clc
close all
rng(1)

addpath('spectral_method_3D')
addpath('chebfun')
addpath('tensor_recursive')
fprintf('Please run fast-poisson-solvers-master/code/setup.m by hand to set the paths for the nested alternating direction implicit method.\n \n')

x = 2*(rand(1000,1)-0.5); y = 2*(rand(1000,1)-0.5); z = 2*(rand(1000,1)-0.5);

%% iterate over different values of n
for nn = [11,31,51,151]
n = nn*[1,1,1];
uTrue = @(x,y,z) sin(pi*x).*sin(pi*y).*sin(pi*z);
f = @(x,y,z) -3*pi*pi*sin(pi*x).*sin(pi*y).*sin(pi*z);
bcFun = @(x,y) 0.*x.*y;

if nn < 40
% run NADIM and meassure the runtime
tic()
F = coeffs3(chebfun3(@(x,y,z)f(x,y,z)),nn,nn,nn);
uNADIM = poisson_cube(F);
timeNADIM = toc();

% error analysis for NADIM
funNADIM = chebfun3( uNADIM, 'coeffs' );
for i = 1:1000
    errNADIM(i) = abs(funNADIM(x(i),y(i),z(i))-uTrue(x(i),y(i),z(i)));
end
errNADIM = max(errNADIM);
fprintf('For n=%3.i NADIM     uses %.2d seconds and achieves an error of %.2d .\n',nn,timeNADIM,errNADIM)
end

% run the recursive solver and meassure the runtime
tic()
LCP{1,1} = [0,0,1]; LCP{1,2} = [1,0,0]; LCP{1,3} = [1,0,0];
LCP{2,1} = [1,0,0]; LCP{2,2} = [0,0,1]; LCP{2,3} = [1,0,0];
LCP{3,1} = [1,0,0]; LCP{3,2} = [1,0,0]; LCP{3,3} = [0,0,1];
[opLCP,lambda] = getForwardOperatorUltra(LCP,n); 
rhsCoeffs = getFullCoeffsFromFunctionUltra(f,n,lambda);
[T1,F1,T2,F2,T3,F3] = getBoundaryConditionMatrices(n,bcFun,bcFun,bcFun,bcFun,bcFun,bcFun);
uRec = solveWithElimination(opLCP,rhsCoeffs,T1,F1,T2,F2,T3,F3,1); %flag for recursive solver
timeRec = toc();

% error analyis for recursive solver
for i = 1:1000
    errRec(i) = abs(funeval(uRec,x(i),y(i),z(i))-uTrue(x(i),y(i),z(i)));
end
errRec = max(errRec);
fprintf('For n=%3.i recursive uses %.2d seconds and achieves an error of %.2d the matrix L1y has condition %.2d.\n',nn,timeRec,errRec,condest(opLCP{1,2}))

if nn < 60
% run backslash after reshape and meassure the runtime
tic()
LCP{1,1} = [0,0,1]; LCP{1,2} = [1,0,0]; LCP{1,3} = [1,0,0];
LCP{2,1} = [1,0,0]; LCP{2,2} = [0,0,1]; LCP{2,3} = [1,0,0];
LCP{3,1} = [1,0,0]; LCP{3,2} = [1,0,0]; LCP{3,3} = [0,0,1];
[opLCP,lambda] = getForwardOperatorUltra(LCP,n); 
rhsCoeffs = getFullCoeffsFromFunctionUltra(f,n,lambda);
[T1,F1,T2,F2,T3,F3] = getBoundaryConditionMatrices(n,bcFun,bcFun,bcFun,bcFun,bcFun,bcFun);
uReshape = solveWithElimination(opLCP,rhsCoeffs,T1,F1,T2,F2,T3,F3,0); %flag for reshape
timeReshape = toc();

% error analyis for backslash after reshape
for i = 1:1000
    errReshape(i) = abs(funeval(uReshape,x(i),y(i),z(i))-uTrue(x(i),y(i),z(i)));
end
errReshape = max(errReshape);
fprintf('For n=%3.i reshape   uses %.2d seconds and achieves an error of %.2d .\n',nn,timeReshape,errReshape)
end
end



