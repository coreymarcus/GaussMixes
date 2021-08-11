function [n] = UnitVectorExp(phi)
%UnitVectorExp

%locals
nPhi = norm(phi);

if(nPhi == 0)
    n = [0 0 1]';
else

%create output
n = [sin(nPhi)/nPhi*phi; cos(nPhi)];
end


end