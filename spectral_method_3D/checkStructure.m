function [LCP,FastFlag] = checkStructure(LCP)
% Tests whether the CP decomposition in LCP leads to a structure, which we
% can exploit in the fast solver. The returned LCP is permuted to have the
% desired symmetries.

FastFlag = 0;

if size(LCP,1) ~= 3 || size(LCP,2) ~= 3
    return
end

% find same pairs in first mode 
diff1 = zeros([3,1]);
a = LCP{1,1};
b = LCP{1,2};
c = LCP{1,3};
if b == c 
    diff1(1) = 1;
end
if a == c
    diff1(2) = 1;
end
if a == b 
    diff1(3) = 1;
end
if sum(diff1) < 1 
    return
end

% find same pairs in second mode 
diff2 = zeros([3,1]);
a = LCP{2,1};
b = LCP{2,2};
c = LCP{2,3};
if b == c 
    diff2(1) = 1;
end
if a == c
    diff2(2) = 1;
end
if a == b 
    diff2(3) = 1;
end
if sum(diff2) < 1 
    return
end

% find same pairs in first mode 
diff3 = zeros([3,1]);
a = LCP{3,1};
b = LCP{3,2};
c = LCP{3,3};
if b == c 
    diff3(1) = 1;
end
if a == c
    diff3(2) = 1;
end
if a == b 
    diff3(3) = 1;
end
if sum(diff3) < 1 
    return
end

if min(diff1 + diff2 + diff3) == 0    
    return
end

% there has to be a permutation allowing the use of the fast solver
P = perms([1,2,3]);


for iter = 1:size(P,1)
    p = P(iter,:);
    if min([diff1(p(1)),diff2(p(2)),diff3(p(3))]) > 0
       % we found the permutation
        LCP2 = LCP;
        LCP2{1,p(1)} = LCP{1,1};
        LCP2{1,p(2)} = LCP{1,2};
        LCP2{1,p(3)} = LCP{1,3};
        LCP2{2,p(1)} = LCP{2,1};
        LCP2{2,p(2)} = LCP{2,2};
        LCP2{2,p(3)} = LCP{2,3};
        LCP2{3,p(1)} = LCP{3,1};
        LCP2{3,p(2)} = LCP{3,2};
        LCP2{3,p(3)} = LCP{3,3};
        LCP = LCP2;
        FastFlag = 1;
        break 
    end
end



end

