function X = solveSylvester(opLCP,D,fastFlag)
% solves the tensor-valued linear equation using either the fast recursive solver for
% Laplace-like equations or by reshaping and using backslash

r = size(opLCP,2);

% check if we can use the fast recursive solver
if r == 3 && fastFlag == 1
    % use the fast recursive solver
    AA = opLCP{1,2}\opLCP{1,1};
    BB = opLCP{2,1}\opLCP{2,2};
    CC = opLCP{2,1}\opLCP{3,3};
    
    DD = invtprod(D,opLCP{1,2},opLCP{2,1},opLCP{3,1});
    X = real(laplace_recursive( {full(AA),full(BB),full(CC)}, DD ));
else
    % use reshape and backslash to solve when the fast solver can not be applied
    A = getForwardOperatorMatrix(opLCP);
    X = reshape(A\D(:),size(D));
end
end