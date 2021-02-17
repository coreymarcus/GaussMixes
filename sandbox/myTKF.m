clear
close all
clc

%parameters
mu_x = 20;
var_x = 4;
mu_z = 0;
var_z = 1;
N = 4000;
mu_w = 0;
var_w = 1;

%draw initial estimate of x
xhatTKF = zeros(1,N);
PhatTKF = zeros(1,N);
xhatTKF(1) = mvnrnd(mu_x, var_x);
PhatTKF(1) = var_x;
xhatKF = zeros(1,N);
PhatKF = zeros(1,N);
xhatKF(1) = xhatTKF(1);
PhatKF(1) = var_x;

%cycle
for ii = 2:N
    
    %draw noise
    z = mvnrnd(mu_z,var_z);
    w = mvnrnd(mu_w,var_w);
    
    %measurement
    y = (1 + z)*mu_x + w;
    
    %variances
    PxyTKF = PhatTKF(ii-1);
    PyyTKF = PhatTKF(ii-1) + var_w + PhatTKF(ii-1)*var_z + var_z*xhatTKF(ii-1)^2;
    PxyKF = PhatKF(ii-1);
    PyyKF = PhatKF(ii-1) + var_w;
    
    %kalman gain
    KTKF = PxyTKF/PyyTKF;
    KKF = PxyKF/PyyKF;
    
    %update
    xhatTKF(ii) = xhatTKF(ii-1) + KTKF*(y - xhatTKF(ii-1));
    PhatTKF(ii) = (1 - KTKF)*PhatTKF(ii-1)*(1 - KTKF) + KTKF*var_w*KTKF;
    xhatKF(ii) = xhatKF(ii-1) + KKF*(y - xhatKF(ii-1));
    PhatKF(ii) = (1 - KKF)*PhatKF(ii-1)*(1 - KKF) + KKF*var_w*KKF;
    
    
end


figure
subplot(2,1,1)
plot(xhatTKF)
hold on
plot(xhatKF)
subplot(2,1,2)
plot(PhatTKF)
hold on
plot(PhatKF)



