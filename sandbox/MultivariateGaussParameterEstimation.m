clear
close all
clc

%% Setup

%seed
rng(4)

n = 10; %number of measurements with each cycle
N_dist = 100000; %number of samples for checking distributions
N_cycle = 100; %number of estimation cycles

%true distribution
p = 2;
mu_true = [3; 2];
var_true = .5*eye(p);

%initial estimate of distribution
P_mu0 = .4*eye(p);
P_var0 = .1*eye(p*p);
mu_hat0 = mvnrnd(mu_true,P_mu0);
var_hat0 = mvnrnd(var_true,P_var0);

%% Check distributions

%initialize storage
mu_meas = zeros(N_dist,1);
var_meas = zeros(N_dist,1);

for ii = 1:N_dist
    
    %grab samples of true distribution
    y = mvnrnd(mu_true, var_true,n);
    
    %get measurements of mu and var
    mu_meas(ii) = mean(y);
    var_meas(ii) = sum((y - mean(y)).^2)/n; %MMSE
%     var_meas(ii) = sum((y - mean(y)).^2)/(n-1); %unbiased
%     var_meas(ii) = sum((y - mu_true).^2)/(n-1); %uncorrelated
    
end

%% Estimation

%initialize storage
mu_hat = zeros(N_cycle + 1,1);
var_hat = zeros(N_cycle + 1, 1);
P_mu = zeros(N_cycle + 1,1);
P_var = zeros(N_cycle + 1, 1);
n_meas = zeros(N_cycle + 1,1);
nu_meas = zeros(N_cycle + 1,1);
prodterm = zeros(N_cycle + 1,1);
mu_hat(1) = mu_hat0;
var_hat(1) = var_hat0;
P_mu(1) = P_mu0;
P_var(1) = P_var0;
n_meas(1) = 1;
nu_meas(1) = 1;
prodterm(1) = nu_meas(1)*var_hat(1);

for ii = 1:N_cycle
    
    %grab samples of true distribution
    y = mvnrnd(mu_true, var_true,n);
    
    % My algorithm (definitely flawed)
    %     %last cycle estimates
    %     mu_bar = mu_hat(ii);
    %     %     var_bar = var_hat(ii);
    %     var_bar = var_true;
    %     Pmu_bar = P_mu(ii);
    %     alpha_bar = alpha_var(ii);
    %     beta_bar = beta_var(ii);
    %     
    %     %update estimate of mu
    %     R_mu = var_bar/n;
    %     mu_hat(ii+1) = mu_bar + (Pmu_bar/(Pmu_bar + R_mu))*(mu_meas(ii) - mu_bar);
    %     P_mu(ii+1) = Pmu_bar - Pmu_bar^2/(Pmu_bar + R_mu);
    %     
    %     %update estimate of var
    %     alpha_meas = (n-1)/2;
    % %     beta_meas = n/(2*var_bar); %MMSE
    %     beta_meas = (n-1)/(2*var_bar); %unbiased
    %     alpha_var(ii+1) = alpha_bar + alpha_meas;
    %     beta_var(ii+1) = beta_bar + beta_meas/var_meas(ii);
    %     var_hat(ii+1) = alpha_var(ii+1)/beta_var(ii+1);
    %     P_var(ii+1) = alpha_var(ii+1)/beta_var(ii+1)^2;
    
    % Wikipedia Algorithm
    % (https://en.wikipedia.org/wiki/Normal_distribution#Estimation_of_parameters)
    
    %last cycle numbers
    xbar = mean(y);
    n0 = n_meas(ii);
    nu0 = nu_meas(ii);
    mu0 = mu_hat(ii);
    prod0 = prodterm(ii);
    
    %update
    mu_hat(ii+1) = (n0*mu0 + n*xbar)/(n0 + n);
    n_meas(ii+1) = n0 + n;
    nu_meas(ii+1) = nu0 + n;
    prodterm(ii+1) = prod0 + sum((y - xbar).^2) + (n0*n)*(mu0 - xbar)^2/(n0 + n);
    
    % human understandable estimates
    P_mu(ii+1) = prodterm(ii+1)/(nu_meas(ii+1)*n_meas(ii+1));
    var_hat(ii+1) = (0.5*prodterm(ii+1))/(0.5*nu_meas(ii+1) - 1);
    P_var(ii+1) = (0.5*prodterm(ii+1))^2/( (0.5*nu_meas(ii+1) - 1)^2 * (0.5*nu_meas(ii+1) - 2));
    
end

%% Plotting

%sample what we expect from the distributions
mu_sample = mvnrnd(mu_true, var_true/n, N_dist);
var_sample = gamrnd((n-1)/2,(2*var_true)/n,N_dist,1); %MMSE
% var_sample = gamrnd((n-1)/2,(2*var_true)/(n-1),N_dist,1); %unbiased

Nbins = 100;

figure
histogram(mu_meas,Nbins,'Normalization','probability','DisplayStyle','stairs')
title('Measurement of Mean')
hold on
histogram(mu_sample,Nbins,'Normalization','probability','DisplayStyle','stairs')
legend('Distribution of Measurements','Expected Distribution of Measurements')

figure
histogram(var_meas,Nbins,'Normalization','probability','DisplayStyle','stairs')
title('Measurement of Var')
hold on
histogram(var_sample,Nbins,'Normalization','probability','DisplayStyle','stairs')
legend('Distribution of Measurements','Expected Distribution of Measurements')


figure
subplot(2,1,1)
plot(mu_hat - mu_true)

subplot(2,1,2)
plot(var_hat - var_true)

figure
subplot(2,1,1)
plot(P_mu)

subplot(2,1,2)
plot(P_var)