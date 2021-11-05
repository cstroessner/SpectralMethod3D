function v = evalChebyPolys(n,x)
% evaluates the first n Chebyshev polynomials at x

for i = 1:n
    v(i) = cos((i-1)*acos(x));
end
v = v';
end
