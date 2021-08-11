function [phi] = UnitVectorSubtract(n,m)
%UnitVectorSubtract

% locals
d = n'*m;
lambda_m = m(3);

if ((d >= 0) && (lambda_m >= 0))
    prod =  UnitVectorToPlus(m)'*n;
elseif ((d >= 0) && (lambda_m < 0))
    prod =  UnitVectorToPlus(-m)'*-n;    
elseif ((d < 0) && (lambda_m >= 0))
    prod =  -1*(UnitVectorToPlus(-m)'*-n);    
elseif ((d < 0) && (lambda_< 0))
    prod =  -UnitVectorToPlus(m)'*n;    
else
    error('Bad case!')    
end

phi = UnitVectorLog(prod);


end

