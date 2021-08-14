function [mPerp] = PlanePerpVec(n, d, m)
%PlanePerpVec

%ensure column
if(size(n,2) > 1)
    n = n';
end

if(size(m,2) > 1)
    m = m';
end

%check length
if(abs(norm(n) - 1) > 0.0001)
    warning('Non-unit vector input')
end

mPerp = (n'*m - d)*n;

end

