function D = CPdec(T)
% Computes a CP decomposition of the tensor T naively.  

if ndims(T) ~= 3
   warning('Dimension of coefficient tensor does not fit.') 
   return
end

% naive approach
n = size(T);
I = 1;
for i = 1:n(1)
    for j = 1:n(2)
        for k = 1:n(3)
            if T(i,j,k) ~= 0
               D{1,I} = (1:n(1) == i)*T(i,j,k);
               D{2,I} = 1:n(2) == j;
               D{3,I} = 1:n(3) == k;
               I = I+1;
            end
        end
    end   
end
    
    
end

