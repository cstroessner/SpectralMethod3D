function [opLCP,lambda] = getForwardOperatorUltraNonConst(Split,n) 
% Computes the forward differential operator as map from 3D Chebyshev
% coefficients to 3D C^{lambda} ultraspherical coefficients. The returned matrices need to
% be multiplied in the corresponding modes and summed up afterwards. See
% applyForwardOperator for details.

lambda = [size(Split{1,1},2)-1,size(Split{2,1},2)-1,size(Split{3,1},2)-1];
r = size(Split,1);
opLCP = cell(size(Split));

for alpha = 1:r
    B = sparse(n(1),n(1));
    for k = 1:lambda(1)+1
        A = getDUltra(n(1),k-1);
        M = MultiplicationMatrix(Split{alpha,1}(:,k), k-1);
        A = M* A;
        for l = k:lambda(1)
            A = getSUltra(n(1),l-1) * A;
        end
        B = B+A;
    end
    opLCP{1,alpha} = B;
end
for alpha = 1:r
    B = sparse(n(2),n(2));
    for k = 1:lambda(2)+1
        A = getDUltra(n(2),k-1);
        M = MultiplicationMatrix(Split{alpha,2}(:,k), k-1);
        A = M* A;
        for l = k:lambda(2)
            A = getSUltra(n(2),l-1) * A;
        end
        B = B+A;
    end
    opLCP{2,alpha} = B;
end
for alpha = 1:r
    B = sparse(n(3),n(3));
    for k = 1:lambda(3)+1
        A = getDUltra(n(3),k-1);
        M = MultiplicationMatrix(Split{alpha,3}(:,k), k-1);
        A = M * A;
        for l = k:lambda(3)
            A = getSUltra(n(3),l-1) * A;
        end
        B = B+A;
    end
    opLCP{3,alpha} = B;
end
end