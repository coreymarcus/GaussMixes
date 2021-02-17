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

%use k-means to find initial guess
[mu_hat, var_hat, w_hat] = kMeans(2, y);

%run standard EM
options.loopmax = 100;
options.tol = 1E-8;
model.mu = mu_hat(:,1);
model.var = var_hat(:,1);
model.w = w_hat(:,1);
model.k = 2;
EMmodel = EM(model, y, options);

%run EM_MAP
MAP_model = EM_MAP(model, y, options);

%draw another set of measurements
N2 = 2;
y2 = zeros(N2,1);
for ii = 1:N2
    %which distribution do we draw from
    idx = randi(2,1);
    
    %draw
    y2(ii) = mvnrnd(mu_true(idx),var_true(idx));
end

%run EM_MAP again
MAP_model2 = EM_MAP(MAP_model, y2, options);