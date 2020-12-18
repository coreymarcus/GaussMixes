%Script to test TLS function
clear
close all
clc

addpath("../matlabScripts")

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

% %detirministic points
% x_meas = [3 4 5 6 7 8 9];
% y_meas = [7 7 11 11 15 16 19];

%TLS
[m_hat_TLS, b_hat_TLS] = TLS(x_meas,y_meas);

%least squares
H = ones(Ninit,2);
H(:,1) = x_meas;
x_hat_LS = inv(H'*H)*H'*y_meas;
m_hat_LS = x_hat_LS(1);
b_hat_LS = x_hat_LS(2);

%sample
x_sample = [x1 x2];
y_sample_TLS = m_hat_TLS*x_sample + b_hat_TLS;
y_sample_LS = m_hat_LS*x_sample + b_hat_LS;
y_true = m_true*x_sample + b_true;

%plot
figure
plot(x_sample,y_true,x_sample,y_sample_TLS,x_sample,y_sample_LS)
hold on
scatter(x_meas,y_meas,'x')
legend('True','TLS','LS','Data','Location','best')
xlabel('x')
ylabel('y')