%Try to estimate the equation of a line

clear
close all
clc

addpath("../matlabScripts")

rng(3)

%setup truth
m_true = 2;
b_true = 1;

%x domain
x1 = 1;
x2 = 4;

%variances of measurement noises
sig2 = .1;

%draw some initial points
Ninit = 10;
x_init = (x2 - x1)*rand(Ninit,1) + x1;
y_init = m_true*x_init + b_true;

%corrupt with noise
x_meas = mvnrnd(x_init, sig2*eye(Ninit))';
y_meas = mvnrnd(y_init, sig2*eye(Ninit))';

%perform an initial total least squares estimate
[m_hat, b_hat, Phat] = TLS(x_meas,y_meas, sig2, sig2);
xhat = [m_hat;
    b_hat];

%plot
Nsample = 10;
x_sample = linspace(x1,x2,Nsample);
y_sample = m_hat*x_sample + b_hat;
figure
scatter(x_meas,y_meas,'rx');
hold on
plot(x_sample,y_sample,'b')

%begin estimation
Ncycle = 100;
Ndraw = 2;
for ii = 1:Ncycle
    
    %draw truth
    x_true = (x2 - x1)*rand(Ndraw,1) + x1;
    y_true = m_true*x_true + b_true;
    
    %corrupt with noise
    x_draw = mvnrnd(x_true, sig2*eye(Ndraw))';
    y_draw = mvnrnd(y_true, sig2*eye(Ndraw))';
    x_meas = [x_meas; x_draw];
    y_meas = [y_meas; y_draw];
    
    %grab old estimates
    xbar = xhat(:,end);
    Pbar = Phat(:,:,end);
    
    %create H
    H = ones(Ndraw,2);
    H(:,1) = x_draw;
    
    %create total dispersion matrix
    Q = blkdiag(sig2*eye(Ndraw),zeros(Ndraw),sig2*eye(Ndraw));
    
    %Call weighted total Kalman Filter
    [xhat_new, Phat_new] = WTKF(y_draw, xbar, Pbar, Q, H);
    
    xhat(:,end+1) = xhat_new;
    Phat(:,:,end+1) = Phat_new;
    
end

%perform an initial total least squares estimate
[m_hat_final, b_hat_final, ~] = TLS(x_meas,y_meas, sig2, sig2);

%plot
Nsample = 10;
x_sample = linspace(x1,x2,Nsample);
y_sampleBatch = m_hat_final*x_sample + b_hat_final;
y_sampleKF = xhat(1,end)*x_sample + xhat(2,end);
y_sampleInit = xhat(1,1)*x_sample + xhat(2,1);
figure
scatter(x_meas,y_meas,'rx');
hold on
plot(x_sample,y_sampleBatch,'b',x_sample,y_sampleKF,'k',...
    x_sample,y_sampleInit,'g')
legend('data','batch','KF','Initial','location','best')

figure
subplot(2,1,1)
plot(xhat(1,:))
title("Estimate")
ylabel('m')
grid on
subplot(2,1,2)
plot(xhat(2,:))
grid on
ylabel('b')

figure
subplot(2,1,1)
semilogy(squeeze(Phat(1,1,:)))
title('Variance')
ylabel('var(m)')
subplot(2,1,2)
semilogy(squeeze(Phat(2,2,:)))
ylabel('var(b)')