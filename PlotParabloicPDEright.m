% this generates the plot in Figure 4b)

addpath('spectral_method_3D')
addpath('tensor_recursive')
clear
clc
close all
rng(1)
tic();
Col = {'b','r','m'};

% initialize params
N = 20;
h = 0.0001;
tsteps= 2000;
n = [N,N,N];

% set Functions
utrue = @(x,y,z,t) exp(-3*pi*pi*t).*sin(pi*x).*sin(y*pi).*sin(z.*pi);
uStart = @(x,y,z) utrue(x,y,z,0);


%% Solve with chebop3

% initialize test points
x = rand([100,1])*2-1;
y = rand([100,1])*2-1;
z = rand([100,1])*2-1;

% discretize the functions
u = getFullCoeffsFromFunction(uStart,n);

% discretize the PDE corresponding to (I-hD) x = u; 
LCP = cell(3);
LCP{1,1} = [1,0,-h];
LCP{1,2} = [1,0,0];
LCP{1,3} = [1,0,0];
LCP{2,1} = [1,0,0];
LCP{2,2} = [0,0,-h];
LCP{2,3} = [1,0,0];
LCP{3,1} = [1,0,0];
LCP{3,2} = [1,0,0];
LCP{3,3} = [0,0,-h];
[LCP,fastFlag] = checkStructure(LCP);
[opLCP,lambda] = getForwardOperatorUltra(LCP,n);
bc = @(x,y) 0.*x.*y;
[T1,F1,T2,F2,T3,F3] = getBoundaryConditionMatrices(n,bc,bc,bc,bc,bc,bc);

% implicit Euler
for step = 1:tsteps
    
    step
    
    tic()
    
    % solve (I-hL) x = u;
    % step 1 set up rhs in terms of ultrasphericals
    rhs = tprod(u,getSUltra(n(1),0),getSUltra(n(2),0),getSUltra(n(3),0));
    rhs = tprod(rhs,getSUltra(n(1),1),getSUltra(n(2),1),getSUltra(n(3),1));
    
    % step 2 solve the PDE including the boundary conditions of u
    u = solveWithElimination(opLCP,rhs,T1,F1,T2,F2,T3,F3,fastFlag);
    
    tChebop3(step) = toc();
    
    % error analysis
    t = step*h;
    for i = 1:100
        trueVal(i) = utrue(x(i),y(i),z(i),t);
        approxVal(i) = funeval(u,x(i),y(i),z(i));
    end
    abserr(step) = max(abs(trueVal-approxVal));
    relerr(step) = abserr(step)/exp(-3*pi*pi*t);
  
end

%% Finite differences for comparison

% discretization params
ht = h;
hx = 2/(N-1);

% discretize function
nodes = linspace(-1,1,N);
[X,Y,Z]=ndgrid(nodes,nodes,nodes);
U = uStart(X,Y,Z);

% assemble differential operator
e = ones(N,1);
tmpD = (1/hx)^2*spdiags([e -2*e e],-1:1,N,N);
D = sparse(N,N);
D(2:N-1,:) = tmpD(2:N-1,:);
I = speye(N,N);
D = kron(D,kron(I,I)) + kron(I,kron(D,I)) + kron(I,kron(I,D));

% solve with implicit Euler
for step = 1:tsteps
    
    step
    
    U = U(:);
    tic()
    U = (speye(N^3)-ht*D)\U;
    tFiniteDiff(step) = toc();
    U = reshape(U,n);
    t = step*ht;
    for i = 1:100
        trueVal(i) = utrue(x(i),y(i),z(i),t);
    end
    approxVal = interp3(nodes,nodes,nodes,U,x,y,z)';
    abserrFD(step) = max(abs(trueVal-approxVal));
    relerrFD(step) = abserrFD(step)/exp(-3*pi*pi*t);
end

%% plot error
set(gca,'fontsize',10)
cols = 1;
set(figure(1), 'Position', [0 0 470 400])
semilogy(1:tsteps,abserrFD,'color',Col{cols}) 
hold on
semilogy(1:tsteps,relerrFD,'color',Col{cols},'linestyle',':')
cols = cols+1;
semilogy(1:tsteps,abserr,'color',Col{cols})
semilogy(1:tsteps,relerr,'color',Col{cols},'linestyle',':')
xlabel('time step $\tau$','Interpreter','latex')
ylabel('estimated error','Interpreter','latex')
cols = cols+1;
leg = legend('absolute error (FD)','relative error (FD)','absolute error (GSM)','relative error (GSM)')
set(leg,'Interpreter','latex','Location','southeast');
%print -depsc 'TimeDependentSmallH'

tChebop3 = sum(tChebop3)
tFiniteDiff = sum(tFiniteDiff)


