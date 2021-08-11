function [nPlus] = UnitVectorToPlus(n)
%UnitVectorToPlus

%ensure column
if(size(n,2) > 1)
    n = n';
end

%locals
kappa = n(1:2);
lambda = n(3);

%check length
% if(abs(norm(n) - 1) > 0.0001)
%     warning('Non-unit vector input')
% end

%poutput
nPlus = [[1 0;0 1] - kappa*kappa'/(1+lambda), kappa;
    -kappa', lambda];
end

