%This script performs some simple Expectation Maximization for learning
%purposes

clear
close all
clc

addpath("../matlabScripts")

%% Initialization

%seed
rng(3)

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

%plot measurements
figure
hist(y,50)

%% E-M

%use k-means to find initial guess
[mu_hat, var_hat, w_hat] = kMeans(2, y);

%initialize gamma's and n's
gamma = zeros(N,2); %prob that (row) measurement from (col) gaussian
n = zeros(2,1); %intermediate variable

loopcounter = 1;
looplogic = true;
loopmax = 100;
while looplogic
    
    % E-step
    
    %cycle through measurements
    for ii = 1:N
        
        %cycle through gaussians
        for jj = 1:2
            gamma(ii,jj) = w_hat(jj,end)*gaussEval(y(ii),mu_hat(jj,end),var_hat(jj,end));
        end
        
        %normalize
        gamma(ii,:) = gamma(ii,:)/sum(gamma(ii,:));
    end
    
    %calculate n
    n(1) = sum(gamma(:,1));
    n(2) = sum(gamma(:,2));
    
    % M-step
    w_hat = [w_hat, n/N];
    mu_hat = [mu_hat, [sum(gamma(:,1).*y)/n(1);
        sum(gamma(:,2).*y)/n(2)]];
    var_hat = [var_hat, [sum(gamma(:,1).*(y - mu_hat(1,end)).^2)/n(1);
        sum(gamma(:,2).*(y - mu_hat(2,end)).^2)/n(2)]];
    
    %loop management
    loopcounter = loopcounter + 1;
    if(loopcounter >= loopmax)
        looplogic = false;
    end
    
end

%run EM function to see if we get the same thing
options.loopmax = 100;
options.tol = 1E-8;
model.mu = mu_hat(:,1);
model.var = var_hat(:,1);
model.w = w_hat(:,1);
model.k = 2;
[EMmodel] = EM(model,y, options);

figure
subplot(3,1,1)
plot(mu_hat(1,:))
hold on
plot(mu_hat(2,:))


subplot(3,1,2)
plot(var_hat(1,:))
hold on
plot(var_hat(2,:))


subplot(3,1,3)
plot(w_hat(1,:))
hold on
plot(w_hat(2,:))
