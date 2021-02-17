%Test of WTKF based on paper

clear
close all
clc

%define constants from paper
a = [1.3, 0.8, 1.0]';
t = 5:5:25;
f = [a*5;
    zeros(3,1)];
x0 = [7.86, 8.74, 14.4, -0.08, 0.28, 0.07]';
n = length(x0);
Sigma0 = zeros(6);
Sigma0(4:6,4:6) = 0.01*eye(3);
Sigma0(1:3,1:3) = [5 5 1;
    5 11 4;
    1 4 2];
theta_i = zeros(6);
theta_i(1:3,1:3) = [2.9, 3.4, 1.0;
    3.4, 6.0, 3.2,;
    1.0, 3.2, 2.4];
theta_i(4:6,4:6) = 0.0025*eye(3);
q = [1 0 0 1 0.1 0 0.4 0.6 0.2;
    0 1 0 0 0 0 0.1 0.2 0.1;
    0 1 1 0.9 0.2 0.7 0.1 0.2 0;
    0.6 0 1 1 1 1 0 0.3 0;
    0.6 0.1 1 0.6 0.1 1 0.6 0.1 1;
    0.2 0.3 0.4 0.6 0.1 1 0.6 0.1 1;
    0.4 0.7 0.8 0.4 0.2 0.4 0.6 0.1 1;
    0.1 0.3 0.5 0.6 0.6 0.4 0.3 0.2 1;
    0.6 0.1 1 0.6 0.1 1 0.6 0.1 1];
Q_y = [3.01 1.32 1.08 2.12 1.67 1.35 0.72 1.02 1.58;
    1.32 1.04 0.13 0.73 0.74 0.40 0.14 0.18 0.72;
    1.08 0.13 3.39 2.86 2.40 2.04 1.98 1.81 2.39;
    2.12 0.73 2.86 4.58 3.15 2.27 1.43 2.23 3.22;
    1.67 0.74 2.40 3.15 3.79 3.05 1.26 1.49 3.77;
    1.35 0.40 2.04 2.27 3.05 3.03 1.16 1.33 2.99;
    0.72 0.14 1.98 1.43 1.26 1.16 1.70 1.17 1.21;
    1.02 0.18 1.81 2.23 1.49 1.33 1.17 1.36 1.50;
    1.58 0.72 2.39 3.22 3.77 2.99 1.21 1.50 3.78];
Phi_i = eye(n); %STM

%measurements
m = 9;
y = zeros(9,5);
y(:,1) = [7.82 11.8 19.3 -8.4 -37 89.9 0.08 -33 121]';
y(:,2) = [14.2 15.2 21.8 6.97 -132 128 24.8 -37.5 157]';
y(:,3) = [25.523 21.69 26.15 41.537 -262.67 172.63 50.492 -41.427 195.99]';
y(:,4) = [28.27 24.36 34.82 90.96 -432 228.6 80.97 -44.3 236.9]';
y(:,5) = [32.707 27.548 35.058 149.51 -650.16 287.02 106 -50.29 269.17]';

%measurement mapping
A = zeros(m,n,5);
A(1:3,1:3,1) = eye(3);
A(4:9,1:6,1) = [2.2897 0.329 -2.1843 0.882 -0.2163 -1.1619;
    -1.0252 -6.4313 -0.680 2.1 1.718 4.2552;
    -4.0252 -0.4312 4.3197 -2.9 0.7185 1.255;
    5.945 -8.325 -0.979 0.420 1.065 3.683;
    0.846 -4.041 -2.874 9.015 -0.460 0.9134;
    -1.9252 1.5687 3.319 -3.9 0.718 1.255];

A(1:3,1:3,2) = eye(3);
A(4:9,1:6,2) = [4.564, 1.40 -3.039, 1.068 5.784, 3.95;
    0.632 -7.478, 1.252, -0.269, 7.46, 6.614;
    -2.36, 4.521, 6.252, -8.269, 4.46, 3.614;
    6.393, -4.690, 0.4162, -1.94, 2.178, 7.131;
    0.23, -0.020, -1.080, 8.1352, 2.049, 2.292;
    -0.267, 5.521, 5.252, -6.269, 3.46, 4.614];

A(1:3,1:3,3) = eye(3);
A(4:9,1:6,3) = [4.584 0.892 -3.85 -2.377 3.020 1.036;
    2.508 -13.73 1.342 1.121 8.760 4.22;
    -0.491 4.269 6.342 -9.878 1.760 1.223;
    8.733 -5.631 0.851 0.0870 1.50 9.682;
    1.482 -1.262 -0.755 8.632 0.220 1.893;
    1.608 4.269 5.342 -4.878 1.760 3.22];

A(1:3,1:3,4) = eye(3);
A(4:9,1:6,4) = [8.118 -1.768 -2.702 0.210 4.332 0.412;
    1.90 -20.20 2.511 -0.7763 13.41 6.022;
    -1.099 3.798 7.511 -14.77 3.417 3.022;
    8.15 -7.33 1.7683 -3.103 1.192 13.73;
    1.668 -3.613 0.0826 6.417 -0.0156 1.740;
    1.000 2.798 6.511 -6.7764 3.417 6.022];    

A(1:3,1:3,5) = eye(3);
A(4:9,1:6,5) = [11.75 2.089 -4.201 -2.72 1.55 1.549;
    1.448 -25.266 1.717 1.175 16.73 3.332;
    -1.551 4.734 6.717 -15.825 3.73 0.332;
    7.754 -7.589 0.599 -0.0959 2.246 13.765;
    2.265 -2.522 -1.1247 7.969 -0.0498 -0.281;
    0.5481 2.734 5.717 -4.825 3.73 4.332];

%initialize storage
x_true = zeros(n,5);
x_hat = zeros(n,5);
x_rhat = zeros(n,5);
Phat = zeros(n,n,5);

%begin algorithm
for ii = 1:5
    
    %propagate truth (also step 4)
    if(ii == 1)
        x_true(:,ii) = Phi_i*x0 + f;
        x_hat(:,ii) = Phi_i*x0 + f;
        Phat(:,:,ii) = Sigma0;
    else
        x_true(:,ii) = Phi_i*x_true(:,ii-1) + f;
        x_hat(:,ii) = Phi_i*x_rhat(:,ii-1) + f;
        Phat(:,:,ii) = Pbar;
    end
    
    %find Q_x, Q
    Mat = kron(eye(6),q);
    Mat(9*ii+1:9*ii+3,:) = 0;
    Mat(:,9*ii+1:9*ii+3) = 0;
    Q_x = Mat*Mat';
    Q = blkdiag(Q_x,Q_y);
    
    %Step 1
    mu_rhat = zeros(n,1);
    M = [eye(m*n, m*n), zeros(m*n, m)];
    N = [zeros(m, n*m), eye(m)];
    Im = eye(m,m);
    A_i = A(:,:,ii);
    theta_i = zeros(n,n); %no process noise since problem is static
    E_rhatAi = zeros(m, n); %I believe this is zero because our predicted peturbation is zero
    
    %Step Two: Iteratively Compute lambda_tilde, E_rhat, mu_rhat
    loopcontinue = true;
    maxloop = 100;
    loopidx = 1;
    tol = 1E-5;
    lamb_tilde = zeros(m,1);
    while loopcontinue
        
        %find lambda
        term1 = N - kron((mu_rhat + x_hat(:,ii))',Im)*M;
        lamb_iter = inv(term1*Q*term1' + A_i*(theta_i + Phi_i*Phat(:,:,ii)*Phi_i')*(A_i - E_rhatAi)')*(y(:,ii) - A_i*x_hat(:,ii));
        E_rhat = Q*term1'*lamb_iter;
        mu_rhat = (theta_i + Phi_i*Phat(:,:,ii)*Phi_i')*(A_i - E_rhatAi)'*lamb_iter;
        
        %update E_rhatAi
        vecE_rhatAi = M*E_rhat;
        E_rhatAi = reshape(vecE_rhatAi,m,n);
        
        
        %loop managment
        loopidx = loopidx + 1;
        if(loopidx > maxloop)
            disp("Warning: max loop achieved")
            loopcontinue = false;
        end
        
        if(norm(lamb_iter - lamb_tilde) <= tol)
            loopcontinue = false;
        end
        
        %update lambda
        lamb_tilde = lamb_iter;
    end
    
    %Step 3
    x_rhat(:,ii) = x_hat(:,ii) + mu_rhat;
    Pbar = kron(theta_i + Phi_i*Phat(:,:,ii)*Phi_i',lamb_tilde')*M*Q*M'*...
        kron(theta_i + Phi_i*Phat(:,:,ii)*Phi_i',lamb_tilde')'...
        + Phi_i*Phat(:,:,ii)*Phi_i';
    
    %Step 4 (done as propagate step above)
end

%plots
figure
subplot(3,1,1)
plot(t,x_true(1,:),t,x_rhat(1,:))
title("Translation")
axis([5 25 10 40])

subplot(3,1,2)
plot(t,x_true(2,:),t,x_rhat(2,:))
axis([5 25 10 30])

subplot(3,1,3)
plot(t,x_true(3,:),t,x_rhat(3,:))
xlabel('Time')
axis([5 25 15 45])

figure
subplot(3,1,1)
plot(t,x_true(4,:),t,x_rhat(4,:))
title("Rotation")
axis([5 25 -1 1])

subplot(3,1,2)
plot(t,x_true(5,:),t,x_rhat(5,:))
axis([5 25 -1 1])

subplot(3,1,3)
plot(t,x_true(6,:),t,x_rhat(6,:))
xlabel('Time')
axis([5 25 -1 1])