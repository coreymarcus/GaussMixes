clear
close all
clc

%random number seed
rng(6);

%path for gauss eval
addpath('C:\Users\cm58349\Github\AdvancedEstimationProject')

%domain
xdomain = [0 1];

%number of gaussians
Ngauss = 25;

%variance in x direction
P_x = .01;

%variance in height
P_h = .01;

%draw estimated gaussian locations along x axis
gauss_x_est = rand(Ngauss,1);

%corrupt nominal gaussian locations
gauss_x_true = mvnrnd(gauss_x_est,P_x*eye(Ngauss))';

%generate true gaussian measurements
gauss_h_true = zeros(Ngauss,1);
for ii = 1:Ngauss
    
    if(gauss_x_true(ii) <= 0.5)
        gauss_h_true(ii) = 0;
    else
        gauss_h_true(ii) = 1;
    end
    
end

%generate measured heights
gauss_h_est = mvnrnd(gauss_h_true,P_h*eye(Ngauss))';

%generate sample points
x_sample = xdomain(1):.01:xdomain(2);
Nsample = length(x_sample);

%true height
h_true = zeros(Nsample,1);
for ii = 1:Nsample
    if(x_sample(ii) > 0.5)
        h_true(ii) = 1;
    end
end

%estimated height
h_est_brute = zeros(Nsample,1);
h_est_opt_init_guess = zeros(Nsample,1);
h_est_opt = zeros(Nsample,1);

%generate sample eval points
h_sample = -1:.001:2;
Neval = length(h_sample);

%covariance matrix
P = zeros(2);
P(1,1) = P_x;
P(2,2) = P_h;

%evaluation vector
gausseval = zeros(Neval,Nsample);

%cycle through each sample
for ii = 1:Nsample
    
    %Brute force estimation of maximum
    %cycle through each height evaluation
    for jj = 1:Neval
        
        %cycle though each gaussian
        for kk = 1:Ngauss
            
            %get mu
            mu = [gauss_x_est(kk); gauss_h_est(kk)];
            
            %get eval point
            evalpoint = [x_sample(ii); h_sample(jj)];
            
            %evaluate
            gausseval(jj,ii) = gausseval(jj,ii) + gaussEval(evalpoint, mu, P);
            
        end
        
    end
    
    %find the maximum of the evaluated gaussians
    [~, I] = max(gausseval(:,ii));
    
    %extract
    h_est_brute(ii) = h_sample(I);
    
    %optimization based estimation of maximum
    
    %get initial guess
    opt_w = zeros(Ngauss,1);
    opt_h = gauss_h_est;
    for jj = 1:Ngauss
        opt_w(jj) = gaussEval([x_sample(ii); gauss_h_est(jj)],...
            [gauss_x_est(jj); gauss_h_est(jj)],P);
    end
    
    %normalize weights
    opt_w = opt_w/sum(opt_w);
    h_est_opt_init_guess(ii) = sum(opt_w.*opt_h);
    
    %all the means
    muArray = [gauss_x_est';
        gauss_h_est'];
    
    %function
    f2min = @(h) -1*gaussSum([x_sample(ii); h], muArray, P);
    %
    %     test = zeros(Neval,1);
    %
    %     for jj = 1:Neval
    %         test(jj) = f2min(h_sample(jj));
    %     end
    
    %perform optimization
%     options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
%         'Display','iter-detailed','FiniteDifferenceStepSize',1E-4);
    %     h_est_opt(ii) = lsqnonlin(f2min, h_est_opt_init_guess(ii),[],[],options);
    options = optimoptions('fmincon','Display','off','OptimalityTolerance',1E-8);
    h_est_opt(ii) = fmincon(f2min, h_est_opt_init_guess(ii),[],[],[],[],[],[],[],options);
    
end

%plot the result
figure
plot(x_sample,h_est_brute,'LineWidth',2)
hold on
plot(x_sample,h_est_opt_init_guess,'LineWidth',2)
plot(x_sample,h_est_opt,'LineWidth',2)
errorbar(gauss_x_est,gauss_h_est,...
    sqrt(P_x)*ones(Ngauss,1),sqrt(P_x)*ones(Ngauss,1),...
    sqrt(P_h)*ones(Ngauss,1),sqrt(P_h)*ones(Ngauss,1),'o')
plot(x_sample,h_true,'LineWidth',2)
xlabel('x')
ylabel('h')
legend('GM Estimated Height','Weighted Average (opt. init. guess)',...
    'Optimization Based Estimate','PC Estimated Height w/ 1 \sigma','True Height')

figure
mesh(x_sample,h_sample,gausseval)
hold on
plot3(x_sample,h_true,1000*ones(Nsample,1),'k','LineWidth',2)
view(2)
xlabel('x')
ylabel('h')

%gaussian sum function
function p = gaussSum(x, mu, P)

%get length of mu
N = size(mu,2);

%initialize p
p = 0;

%cycle
for ii = 1:N
    p = p + gaussEval(x, mu(:,ii), P);
end


end