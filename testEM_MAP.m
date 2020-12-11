%Tests the EM_MAP function
clear
close all
clc

%initialize the true system
mu_true = [1, 2.5];
var_true = [.1, .2];
w_true = [0.5, 0.5];

%draw measurements
N = 1000;
y = zeros(N,1);
for ii = 1:N
    %which distribution do we draw from
    idx = randi(2,1);
    
    %draw
    y(ii) = mvnrnd(mu_true(idx),var_true(idx));
end