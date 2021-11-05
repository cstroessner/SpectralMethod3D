function out = funeval(Coeffs, x, y, z)
% evaluates a 3D Chebyshev polynomial with coefficient tensor U at (x,y,z)
    m = [size(Coeffs)];
    u = evalChebyPolys(m(1),x);
    v = evalChebyPolys(m(2),y);
    w = evalChebyPolys(m(3),z);
    
    out = tprod(Coeffs,u',v',w');
end

function v = evalChebyPolys(n,x)
for i = 1:n
    v(i) = cos((i-1)*acos(x));
end
v = v';
end