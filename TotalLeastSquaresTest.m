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
R = sig2*eye(2);

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
[m_hat_TLS, b_hat_TLS, P_LS] = TLS(x_meas,y_meas,sig2,sig2);

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

%for all points, calculate t, the perpendicular distance to the line and s,
%the distance along the line
t_meas = zeros(Ninit,1);
t_true = zeros(Ninit,1);
s_meas = zeros(Ninit,1);
for ii = 1:Ninit
    
    %find t
    t_meas(ii) = (b_hat_TLS + m_hat_TLS*x_meas(ii) - y_meas(ii))/sqrt(1 + m_hat_TLS^2);
    t_true(ii) = (b_hat_TLS + m_hat_TLS*x_init(ii) - y_init(ii))/sqrt(1 + m_hat_TLS^2);
    
    %find s, this should be simplified somehow
    gamma = y_meas(ii) + x_meas(ii)/m_hat_TLS;
    A = [1, -m_hat_TLS;
        1, 1/m_hat_TLS];
    inter = A\[b_hat_TLS; gamma];
    xinter = inter(2);
    s_meas(ii) = xinter*sqrt(1+m_hat_TLS^2);
    
end

%find midpoint of the line
min_s = min(s_meas);
max_s = max(s_meas);
min_x = min_s/sqrt(1+m_hat_TLS^2);
max_x = max_s/sqrt(1+m_hat_TLS^2);
avg_x = 0.5*(min_x + max_x);
avg_y = avg_x*m_hat_TLS + b_hat_TLS;

%approximate the variance in the t direction
pt_pm = avg_x/sqrt(1+m_hat_TLS^2) - (b_hat_TLS + m_hat_TLS*avg_x - avg_y)*m_hat_TLS*(1+m_hat_TLS^2)^(-1.5);
pt_pb = 1/sqrt(1+m_hat_TLS^2);
var_t = [pt_pm, pt_pb]*P_LS*[pt_pm; pt_pb];

%approximate the variance in the s direction
k_unif = 0.75; %hueristic for scaling uniform variance
var_s = k_unif*(1/12)*(max_s - min_s)^2;

%now, rotate covariance into the cartesian frame
P_st = [var_s, 0;
    0, var_t];
theta = atan(m_hat_TLS);
R_st2xy = [cos(theta), -sin(theta);
    sin(theta), cos(theta)];
P_xy = R_st2xy*P_st*R_st2xy';

%plot
figure
plot(x_sample,y_true,x_sample,y_sample_TLS,x_sample,y_sample_LS)
hold on
scatter(x_meas,y_meas,'x')
scatter(avg_x,avg_y)
plot_elipse(gcf, P_xy, [avg_x; avg_y], 3, 'k', false)
plot_elipse(gcf, P_xy + R, [avg_x; avg_y], 3, 'g', false)
legend('True','TLS','LS','Measurements','Midpoint',...
    '3 \sigma covariance ellipse for Gauss',...
    '3 \sigma covariance ellipse for measurements','Location','best')
xlabel('x')
ylabel('y')