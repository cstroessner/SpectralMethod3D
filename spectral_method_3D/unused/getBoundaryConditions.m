function [ABC,rhsBC] = getBoundaryConditions(n,rightbc,leftbc,topbc,bottombc,frontbc,backbc)
% returns a matrix ABC such that ABC*uCoeffs(:) = rhsBC when the boundary
% conditions are fullfilled

% right
coeffs2D = getFull2DCoeffsFromFunction(rightbc,[n(2),n(3)]);
Arbc = []; frbc = [];
for j = 1:n(2)
    for k = 1:n(3)
        Arbc = [Arbc;kron((1:n(3)==k)',kron((1:n(2)==j)',evalChebyPolys(n(1),1)))'];
        frbc = [frbc;coeffs2D(j,k)];
    end
end
% left
coeffs2D = getFull2DCoeffsFromFunction(leftbc,[n(2),n(3)]);
Albc = []; flbc = [];
for j = 1:n(2)
    for k = 1:n(3)
        Albc = [Albc;kron((1:n(3)==k)',kron((1:n(2)==j)',evalChebyPolys(n(1),-1)))'];
        flbc = [flbc;coeffs2D(j,k)];
    end
end
% top
coeffs2D = getFull2DCoeffsFromFunction(topbc,[n(1),n(3)]);
Atbc = []; ftbc = [];
for i = 1:n(1)
    for k = 1:n(3)
        Atbc = [Atbc;kron((1:n(3)==k)',kron(evalChebyPolys(n(2),1),(1:n(1)==i)'))'];
        ftbc = [ftbc;coeffs2D(i,k)];
    end
end
% bottom
coeffs2D = getFull2DCoeffsFromFunction(bottombc,[n(1),n(3)]);
Abbc = []; fbbc = [];
for i = 1:n(1)
    for k = 1:n(3)
        Abbc = [Abbc;kron((1:n(3)==k)',kron(evalChebyPolys(n(2),-1),(1:n(1)==i)'))'];
        fbbc = [fbbc;coeffs2D(i,k)];
    end
end
% front
coeffs2D = getFull2DCoeffsFromFunction(frontbc,[n(1),n(2)]);
Afbc = []; ffbc = [];
for i = 1:n(1)
    for j = 1:n(2)
        Afbc = [Afbc;kron(evalChebyPolys(n(3),-1),kron((1:n(2)==j)',(1:n(1)==i)'))'];
        ffbc = [ffbc;coeffs2D(i,j)];
    end
end
%back
coeffs2D = getFull2DCoeffsFromFunction(backbc,[n(1),n(2)]);
Abackbc = []; fbackbc = [];
for i = 1:n(1)
    for j = 1:n(2)
        Abackbc = [Abackbc;kron(evalChebyPolys(n(3),1),kron((1:n(2)==j)',(1:n(1)==i)'))'];
        fbackbc = [fbackbc;coeffs2D(i,j)];
    end
end
ABC = {Arbc,Albc,Atbc,Abbc,Afbc,Abackbc};
rhsBC = {frbc,flbc,ftbc,fbbc,ffbc,fbackbc};
end
