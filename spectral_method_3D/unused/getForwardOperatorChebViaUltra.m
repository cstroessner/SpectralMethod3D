function [opL1,opL2,opL3] = getForwardOperatorChebViaUltra(L1,L2,L3,n, lambda)
% Computes the forward differential operator as map from 3D Chebyshev
% coefficients to 3D Chebyshev coefficients. The returned matrices need to
% be multiplied in the corresponding modes and summed up afterwards. See
% applyForwardOperator for details.
% Remarks: at the moment we only accept CP decomposition with exactly rank 3
% This function differs from getForwardOperator, since we here use
% ultraspherical coefficients internally.
for alpha = 1:3
    S{alpha} = eye(n(alpha));
    for l = lambda-1:-1:0
        S{alpha} = inv(getSUltra(n(alpha),l))*S{alpha};
    end
end

for alpha = 1:3
    opL1{alpha} = zeros(n(alpha));
    for k = 1:3
        A = getDUltra(n(alpha),k-1);
        for l = k:lambda
            A = getSUltra(n(alpha),l-1) * A;
        end
        A = L1{alpha}(k).* A;
        opL1{alpha} = opL1{alpha}+S{alpha}*A;
    end
end
for alpha = 1:3
    opL2{alpha} = zeros(n(alpha));
    for k = 1:3
        A = getDUltra(n(alpha),k-1);
        for l = k:lambda
            A = getSUltra(n(alpha),l-1) * A;
        end
        A = L2{alpha}(k).* A;
        opL2{alpha} = opL2{alpha}+S{alpha}*A;
    end
end
for alpha = 1:3
    opL3{alpha} = zeros(n(alpha));
    for k = 1:3
        A = getDUltra(n(alpha),k-1);
        for l = k:lambda
            A = getSUltra(n(alpha),l-1) * A;
        end
        A = L3{alpha}(k).* A;
        opL3{alpha} = opL3{alpha}+S{alpha}*A;
    end
end
end