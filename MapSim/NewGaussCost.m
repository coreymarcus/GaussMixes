function [J] = NewGaussCost(x,z,xhat_old,Phat_old, R, Nmeas_old)
%NewGaussCost finds the cost for choosing a new gauss element

%locals
xhat_new = x(1:2);
D = zeros(2);
D(1,1) = x(3);
D(2,2) = x(4);
v1 = [1, x(5)]';
v1 = v1/norm(v1);
v2 = [-x(5), 1]';
v2 = v2/norm(v2);
V = [v1, v2];
Phat_new = V*D/V;
m = size(z,2);

%evaluate liklihood of measurements
J_z = 0;
for ii = 1:m
    J_z = J_z + 1/gaussEval(z(:,ii),xhat_new,Phat_new + R);    
end

%evaluate cost of change in xhat and Phat
xhat_diff = xhat_new - xhat_old;
Phat_diff = reshape(Phat_new - Phat_old, 4, 1);
J_gauss = 1*Nmeas_old * ( xhat_diff'*xhat_diff + Phat_diff'*Phat_diff);

%total cost
J = J_z + J_gauss;

if(~isreal(J))
    disp('Imag')
end

end

