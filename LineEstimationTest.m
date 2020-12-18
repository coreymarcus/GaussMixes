%Try to estimate the equation of a line

clear
close all
clc

rng(3)

%setup truth
m_true = 2;
b_true = 1;

%x domain
x1 = 1;
x2 = 4;

%variances of measurement noises
sig2 = .1;
R = sig2*eye(4);

%draw some initial points
Ninit = 100;
x_init = (x2 - x1)*rand(Ninit,1) + x1;
y_init = m_true*x_init + b_true;

%corrupt with noise
x_meas = mvnrnd(x_init, sig2*eye(Ninit))';
y_meas = mvnrnd(y_init, sig2*eye(Ninit))';

%we'll artificially inflate the line estimate variance with some fake
%process noise
Q = 0.0001*sig2*eye(2);

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
Ndraw = 200;
for ii = 1:Ncycle
    
    %draw truth
    x_draw = (x2 - x1)*rand(Ndraw,1) + x1;
    y_draw = m_true*x_draw + b_true;
    
    %corrupt with noise
    x_draw = mvnrnd(x_draw, sig2*eye(Ndraw))';
    y_draw = mvnrnd(y_draw, sig2*eye(Ndraw))';
    x_meas = [x_meas; x_draw];
    y_meas = [y_meas; y_draw];
    
    %grab old estimates
    xbar = xhat(:,end);
    mbar = xbar(1);
    bbar = xbar(2);
    Pbar = Phat(:,:,end)+Q;
    
    %create H, z and zbar
    H = zeros(2*Ndraw,2);
    z = zeros(2*Ndraw,1);
    zbar = zeros(2*Ndraw,1);
    for jj = 1:Ndraw
        H(2*jj - 1,1) = -(y_draw(jj) - bbar)/mbar^2;
        H(2*jj - 1,2) = -1/mbar;
        H(2*jj, 1) = x_draw(jj);
        H(2*jj, 2) = 1;
        z(2*jj - 1) = x_draw(jj);
        z(2*jj) = y_draw(jj);
        zbar(2*jj - 1) = (y_draw(jj) - bbar)/mbar;
        zbar(2*jj) = x_draw(jj)*mbar + bbar;
    end
    
    %calculate covariances
    Pxz = H*Pbar;
    Pzz = H*Pbar*H' + sig2*eye(2*Ndraw);
    
    %update
    xhat_new = xbar + Pxz'/Pzz*(z - zbar);
    Phat_new = Pbar - Pxz'/Pzz*Pxz;
    Phat_new = 0.5*(Phat_new + Phat_new');
    
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
figure
scatter(x_meas,y_meas,'rx');
hold on
plot(x_sample,y_sampleBatch,'b',x_sample,y_sampleKF,'k')
legend('data','batch','KF','location','best')

figure
subplot(2,1,1)
plot(xhat(1,:))
grid on
subplot(2,1,2)
plot(xhat(2,:))
grid on

figure
subplot(2,1,1)
semilogy(squeeze(Phat(1,1,:)))
subplot(2,1,2)
semilogy(squeeze(Phat(2,2,:)))