%Try to estimate the equation of a line

clear
close all
clc

addpath("../matlabScripts")

% rng(3)

%setup truth
m_avg = 1;
b_avg = 2;
z_avg = [m_avg; b_avg];
P_z = .1*eye(2);

%draw true values
Ncycle = 10000;
z_true = mvnrnd([m_avg; b_avg],P_z,Ncycle)';

%x domain
x1 = -4;
x2 = 40;

%variances of measurement noises
P_x = .1;
P_y = .1;
R = [P_x, 0;
    0, P_y];

%begin estimation
Ndraw = 100;
z_hat = zeros(2,Ncycle);
P_zhat = zeros(2,2,Ncycle);
e_zhat = zeros(2,Ncycle);
res = zeros(Ndraw,Ncycle);
update = zeros(2,Ncycle);
for ii = 1:Ncycle
    
    %draw true x and y values
    x_true = (x2 - x1)*rand(Ndraw,1) + x1;
    y_true = z_true(1,ii)*x_true + z_true(2,ii);
    
    %corrupt with noise
    meas = zeros(Ndraw,2);
    for jj = 1:Ndraw
        meas(jj,:) = mvnrnd([x_true(jj), y_true(jj)],R);
    end
    
    %create H
    H = ones(Ndraw,2);
    H(:,1) = meas(:,1);
    
    %approximate measurement variance
    R2 = ((P_z(1,1) + m_avg^2)*P_x + P_y) * eye(Ndraw);
    
    %calculate covariances
    Pzx = P_z*H';
    Pzz = H*P_z*H' + R2;
    
    %update
    z_hat(:,ii) = z_avg + Pzx/Pzz*(meas(:,2) - H*z_avg);
    P_zhat(:,:,ii) = P_z - Pzx/Pzz*Pzx';
    
    %error
    e_zhat(:,ii) = z_hat(:,ii) - z_true(:,ii);
    
    %residual
    res(:,ii) = meas(:,2) - H*z_avg;
    
    %update
    update(:,ii) = Pzx/Pzz*(meas(:,2) - H*z_avg);
    
end

%error covariance
P_e = cov(e_zhat');

%plot
figure
subplot(2,1,1)
plot(e_zhat(1,:))
hold on
plot(ones(1,Ncycle)*P_e(1,1))
plot(ones(1,Ncycle)*P_zhat(1,1,1))
legend('Error','Error Variance','Predicted Estimate Variance')
ylabel('e_m')

subplot(2,1,2)
plot(e_zhat(2,:))
hold on
plot(ones(1,Ncycle)*P_e(2,2))
plot(ones(1,Ncycle)*P_zhat(2,2,1))
legend('Error','Error Variance','Predicted Estimate Variance')
ylabel('e_b')

%disp some stats
disp("Average Error:")
disp(mean(e_zhat,2))

disp("Error Variance")
disp(P_e)

disp("Estimate Variance")
disp(P_zhat(:,:,1))

mean(update,2)