clear
close all
clc

%% Setup

%seed
rng(2)

n_draw = 10; %number of measurements with each cycle
N_cycle = 1000; %number of estimation cycles

%true distribution
p = 2;
mu_true = [3; 2];
Sig_true = .5*eye(p);
Sig_true(1,2) = .25;
Sig_true(2,1) = .25;

%measurement noise
R = .1*Sig_true;
% R = zeros(p);

%% Estimation

%draw some samples to get the initial estimate of the distribution
m0 = 3;
y = mvnrnd(mu_true, Sig_true + R, m0)';
mu_hat0 = mean(y,2);
Sig_hat0 = -1*R;
for ii = 1:m0
    Sig_hat0 = Sig_hat0 + (1/m0)*(y(:,ii) - mu_hat0)*(y(:,ii) - mu_hat0)';
end
n0 = 4;
Psi0 = (n0-p-1)*Sig_hat0;

%initialize storage
mu_hat = zeros(p,N_cycle + 1);
Sig_hat = zeros(p,p,N_cycle + 1);
m = zeros(1,N_cycle + 1);
Psi = zeros(p,p,N_cycle + 1);
n = zeros(1, N_cycle + 1);

mu_hat(:,1) = mu_hat0;
Sig_hat(:,:,1) = Sig_hat0;
m(1) = m0;
n(1) = n0;
Psi(:,:,1) = Psi0;

for ii = 1:N_cycle
    
    %grab samples of true distribution
    y = mvnrnd(mu_true, Sig_true + R, n_draw)';
    
    % Wikipedia Algorithm
    % (https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Statistical_inference)
    
    %helpers
    xbar = mean(y,2);
    S = -R;
    for jj = 1:n_draw
        S = S + (1/n_draw)*(y(:,jj) - xbar)*(y(:,jj) - xbar)';
    end
    
    %update mean
    mu_hat(:,ii+1) = (1/(n_draw + m(ii)))*(n_draw*xbar + m(ii)*mu_hat(:,ii));
    m(ii+1) = n_draw + m(ii);
    
    %update covariance
    Psi(:,:,ii+1) = Psi(:,:,ii) + n_draw*S + (n_draw*m(ii)/(n_draw + m(ii)))*(xbar - mu_hat(:,ii))*(xbar - mu_hat(:,ii))';
    n(ii+1) = n(ii) + n_draw;
    
    %new estimate of covariance
    Sig_hat(:,:,ii+1) = (1/(n(ii+1) - p - 1))*Psi(:,:,ii+1);
    
end

%% Plotting

figure
plot(mu_hat(1,:) - mu_true(1))
hold on
plot(mu_hat(2,:) - mu_true(2))
grid on

figure
subplot(2,2,1)
plot(squeeze(Sig_hat(1,1,:) - Sig_true(1,1)))
grid on

subplot(2,2,2)
plot(squeeze(Sig_hat(1,2,:) - Sig_true(1,2)))
grid on

subplot(2,2,3)
plot(squeeze(Sig_hat(2,1,:) - Sig_true(2,1)))
grid on

subplot(2,2,4)
plot(squeeze(Sig_hat(2,2,:) - Sig_true(2,2)))
grid on


