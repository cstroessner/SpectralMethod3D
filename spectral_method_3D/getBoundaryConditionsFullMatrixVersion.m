function [ABC,rhsBC] = getBoundaryConditionsFullMatrixVersion(n,rightbc,leftbc,topbc,bottombc,frontbc,backbc)
% returns matrix ABC and vector rhsBS such that a coefficient tensor U
% satisfying the boundary conditions satisfies matBC*U(:)-rhsBC = 0

cheb = @(i,n) cos((i-1).*pi/(n-1));
xx = cheb(1:n(1),n(1));
yy = cheb(1:n(2),n(3));
zz = cheb(1:n(2),n(3));
% right 
Arbc = []; frbc = [];
for j = 1:n(2)
    for k = 1:n(3)          
        Arbc = [Arbc;kron(evalChebyPolys(n(3),zz(k)),kron(evalChebyPolys(n(2),yy(j)),evalChebyPolys(n(1),1)))'];
        frbc = [frbc;rightbc(yy(j),zz(k))];
    end
end
% left 
Albc = []; flbc = [];
for j = 1:n(2)
    for k = 1:n(3)          
        Albc = [Albc;kron(evalChebyPolys(n(3),zz(k)),kron(evalChebyPolys(n(2),yy(j)),evalChebyPolys(n(1),-1)))'];
        flbc = [flbc;leftbc(yy(j),zz(k))];
    end
end
% top
Atbc = []; ftbc = [];
for i = 1:n(1)
    for k = 1:n(3)          
        Atbc = [Atbc;kron(evalChebyPolys(n(3),zz(k)),kron(evalChebyPolys(n(2),1),evalChebyPolys(n(1),xx(i))))'];
        ftbc = [ftbc;topbc(xx(i),zz(k))];
    end
end
% bottom 
Abbc = []; fbbc = [];
for i = 1:n(1)
    for k = 1:n(3)          
        Abbc = [Abbc;kron(evalChebyPolys(n(3),zz(k)),kron(evalChebyPolys(n(2),-1),evalChebyPolys(n(1),xx(i))))'];
        fbbc = [fbbc;bottombc(xx(i),zz(k))];
    end
end
% front
Afbc = []; ffbc = [];
for i = 1:n(1)
    for j = 1:n(2)          
        Afbc = [Afbc;kron(evalChebyPolys(n(3),-1),kron(evalChebyPolys(n(2),yy(j)),evalChebyPolys(n(1),xx(i))))'];
        ffbc = [ffbc;frontbc(xx(i),yy(j))];
    end
end
%back
Abackbc = []; fbackbc = [];
for i = 1:n(1)
    for j = 1:n(2)          
        Abackbc = [Abackbc;kron(evalChebyPolys(n(3),1),kron(evalChebyPolys(n(2),yy(j)),evalChebyPolys(n(1),xx(i))))'];
        fbackbc = [fbackbc;backbc(xx(i),yy(j))];
    end
end
ABC = [Arbc;Albc;Atbc;Abbc;Afbc;Abackbc];
rhsBC = [frbc;flbc;ftbc;fbbc;ffbc;fbackbc];
end

function v = evalChebyPolys(n,x)
for i = 1:n
    v(i) = cos((i-1)*acos(x));
end
v = v';
end