% Code to generate the plots in Figure 2
addpath('spectral_method_3D')
addpath('tensor_recursive')
clear
close all
rng(1)

% get solution, rhs and boundary functions
a = 5; b = 3; c = 5;
k = @(x) a-b*cos(pi*c*x/2);
utrue = @(x,y,z) exp(-k(x)/c).*cos(pi*a*y/2).*cos(pi*b*z/2);
f = @(x,y,z) (pi^2/4*b^2*(sin(pi*c*x/2)).^2.*exp(-k(x)/c) - pi^2/4*b*c*cos(pi*c*x/2).*exp(-k(x)/c)).*cos(pi*a*y/2).*cos(pi*b*z/2) + ...
      exp(-k(x)/c).*cos(pi*a*y/2).*(-pi^2/4*b^2).*cos(pi*b*z/2) + ...
      exp(-k(x)/c).*(-pi^2/4*a^2).*cos(pi*a*y/2).*cos(pi*b*z/2) +...
      k(x).^2.*(exp(-k(x)/c).*cos(pi*a*y/2).*cos(pi*b*z/2));
dirichletData = @(x,y,z)  utrue(x,y,z);
rightbc = @(y,z) dirichletData(1,y,z);
leftbc = @(y,z) dirichletData(-1,y,z);
topbc = @(x,z) dirichletData(x,1,z);
bottombc = @(x,z) dirichletData(x,-1,z);
frontbc = @(x,y) dirichletData(x,y,-1);
backbc = @(x,y) dirichletData(x,y,1);

% build the linear operator and solve the PDE
Niters = 1; Ns = 11:10:151;
for N = Ns
N
n = [N,N,N];
% get Laplacian first
LCP = cell(3);
LCP{1,1} = [0,0,1];
LCP{1,2} = [1,0,0];
LCP{1,3} = [1,0,0];
LCP{2,1} = [1,0,0];
LCP{2,2} = [0,0,1];
LCP{2,3} = [1,0,0];
LCP{3,1} = [1,0,0];
LCP{3,2} = [1,0,0];
LCP{3,3} = [0,0,1];
[opLCP,lambda] = getForwardOperatorUltra(LCP,n);
% add the kappa^2 u term by modifying opLCP{1,1}
k2Coeffs = get1DCoeffsFromFunction(@(x) k(x).^2,n(1));
opLCP{1,1} = opLCP{1,1} + MultiplicationMatrix(k2Coeffs,2)*getSUltra(n(1),1)*getSUltra(n(1),0);
fastFlag = 1;
% boundary conditions
[T1,F1,T2,F2,T3,F3] = getBoundaryConditionMatrices(n,rightbc,leftbc,topbc,bottombc,frontbc,backbc);
% right hand side
V = getFullCoeffsFromFunctionUltra(f,n,lambda);

% compute solution 
U = solveWithElimination(opLCP,V,T1,F1,T2,F2,T3,F3,fastFlag);

% analyze error
% residual
res(Niters) = computeResidual(U,opLCP,V,T1,F1,T2,F2,T3,F3);
% difference in coefficient tensors
Utrue = getFullCoeffsFromFunction(utrue,n);
errCoeffs(Niters) = max(abs(U-Utrue),[],'all');
% difference in function evaluations
for samples = 1:1000
    x = (2*rand(1)-1); y = (2*rand(1)-1); z = (2*rand(1)-1);
    errFunctionTmp(samples) = abs(utrue(x,y,z)-funeval(U,x,y,z));
    errInterpolationTmp(samples) = abs(utrue(x,y,z)-funeval(Utrue,x,y,z));

end
errEvaluation(Niters) = max(errFunctionTmp);
errInterpolation(Niters) = max(errInterpolationTmp);
conditions(Niters) = condest(opLCP{1,2});
Niters = Niters + 1;
end

%% plot errors
set(gca,'fontsize',10)
cols = 1;
set(figure(1), 'Position', [0 0 470 400])
semilogy(10:10:150,errEvaluation,'r',10:10:150,errInterpolation,'b')
xlabel('n')
ylabel('error')
cols = cols+1;
leg = legend('$||u-u^*||_\infty$','$||\tilde{u}-u^*||_\infty$');
set(leg,'Interpreter','latex','Location','northeast');
xlim([10,150]); ylim([10^(-15),10^0]);
%print -depsc 'HelmholtzError'

% plot solution
figure(2)
[X,Y] = meshgrid(-1:0.05:1,-1:0.05:1);
Z = utrue(X,Y,1/4);
s = surf(X,Y,Z,'FaceAlpha',0.5)
v = [-5 -3 5];
[caz,cel] = view(v)
%print -depsc 'HelmholtzSolution'


