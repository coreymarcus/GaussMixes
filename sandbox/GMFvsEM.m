%I'm going to try and compare EM to a GMF based mapping

clear
close all
% clc

%random number seed
rng(7);

%path for gauss eval
addpath('..\matlabScripts')

%domain
xdomain = [0 1];

%number of gaussians
Ngauss = 2;

%true step mag
step_mag = 0.5;

%number of initially drawn points
Ninit = 100;

%measurement variance in x direction
P_x = .01;

%measurement variance in height
P_h = .0001;

%measurement variance
R = blkdiag(P_x,P_h);

%Number of measurements drawn on each cycle
Ndraw = 50;

%we'll artificially inflate the line estimate variance with some fake
%process noise
Q = 1.0*R;

%number of update cycles
Ncycle = 10;

%draw the true initial point locations
x_init_true = rand(Ninit,1);

%generate true height values
h_init_true = zeros(Ninit,1);
% pop1true = [];
% pop2true = [];
for ii = 1:Ninit
    
    if(x_init_true(ii) <= 0.5)
        h_init_true(ii) = 0;
        %         pop1true = [pop1true;
        %             x_init_true(ii), h_init_true(ii)];
    else
        h_init_true(ii) = step_mag;
        %         pop2true = [pop2true;
        %             x_init_true(ii), h_init_true(ii)];
    end
    
end

%generate measured x and height values
x_init_meas = mvnrnd(x_init_true,P_x*eye(Ninit))';
h_init_meas = mvnrnd(h_init_true,P_h*eye(Ninit))';

%use kMeans to find the initial estimate
data = [x_init_meas';
    h_init_meas'];
[mu_hat_init, P_hat_init, w_hat_init, group_idx] = kMeans(Ngauss, data);

%form EM system
EMmodel.k = Ngauss;
EMmodel.mu = mu_hat_init;
EMmodel.P = P_hat_init;
EMmodel.w = w_hat_init;
EMoptions.loopmax = 100;
EMoptions.tol = 1E-4;

%create a model for the filter
GMmodel.Nmeas = w_hat_init*Ninit; %number of measurements for each mixand
GMmodel.Ntotal = Ninit; %total number of measurements for all mixands
GMmodel.P = zeros(2,2,Ngauss); %mixand covariance
GMmodel.mu = zeros(2,Ngauss); %mixand mean
GMmodel.line = zeros(2,Ngauss); %mixand line slope and intercept
GMmodel.lineP = zeros(2,2,Ngauss); %mixand line slope and intercept covariance
GMmodel.lineMin = zeros(1,Ngauss); %mixand line smallest s value
GMmodel.lineMax = zeros(1,Ngauss); %mixand line largest s value
GMmodel.k = Ngauss; %number of mixands

%use TLS to initialize the estimates of each gaussian
for ii = 1:Ngauss
    xgroup = x_init_meas(group_idx == ii);
    hgroup = h_init_meas(group_idx == ii);
    
    %call TLS
    [m_TLS, b_TLS, P_LS] = TLS(xgroup, hgroup, P_x, P_h);
    GMmodel.line(1,ii) = m_TLS;
    GMmodel.line(2,ii) = b_TLS;
    GMmodel.lineP(:,:,ii) = P_LS;
    
    %we need to project all the points to the closest intersection of the
    %line to find our "first" and "last" points
    s = zeros(length(xgroup),1);
    for jj = 1:length(xgroup)
        gamma = hgroup(jj) + xgroup(jj)/m_TLS;
        A = [1, -m_TLS;
            1, 1/m_TLS];
        inter = A\[b_TLS; gamma];
        xinter = inter(2); %yes, it is the second element
        s(jj) = xinter*sqrt(1+m_TLS^2);
    end
    
    %find min and max s and x values
    min_s = min(s);
    max_s = max(s);
    min_x = min_s/sqrt(1+m_TLS^2);
    max_x = max_s/sqrt(1+m_TLS^2);
    min_y = min_x*m_TLS + b_TLS;
    max_y = max_x*m_TLS + b_TLS;
    avg_x = 0.5*(min_x + max_x);
    avg_y = 0.5*(min_y + max_y);
    GMmodel.lineMin(ii) = min_s;
    GMmodel.lineMax(ii) = max_s;
    GMmodel.mu(:,ii) = [avg_x; avg_y];
    
    %approximate the variance in the t direction
    pt_pm = avg_x/sqrt(1+m_TLS^2) - (b_TLS + m_TLS*avg_x - avg_y)*m_TLS*(1+m_TLS^2)^(-1.5);
    pt_pb = 1/sqrt(1+m_TLS^2);
    var_t = [pt_pm, pt_pb]*P_LS*[pt_pm; pt_pb];
    
    %approximate the variance in the s direction
    k_unif = 0.75; %hueristic for scaling uniform variance
    var_s = k_unif*(1/12)*(max_s - min_s)^2;
    
    %now, rotate covariance into the cartesian frame
    P_st = [var_s, 0;
        0, var_t];
    theta = atan(m_TLS);
    R_st2xy = [cos(theta), -sin(theta);
        sin(theta), cos(theta)];
    GMmodel.P(:,:,ii) = R_st2xy*P_st*R_st2xy';
end

%plot the initial estimate
figure
scatter(x_init_meas(group_idx == 1),h_init_meas(group_idx == 1),'xr')
hold on
scatter(x_init_meas(group_idx == 2),h_init_meas(group_idx == 2),'xb')
for ii = 1:Ngauss
    plot_elipse(gcf, GMmodel.P(:,:,ii), GMmodel.mu(:,ii),...
        GMmodel.Nmeas(ii)*GMmodel.k/GMmodel.Ntotal, '', false)
    plot_elipse(gcf, EMmodel.P(:,:,ii), EMmodel.mu(:,ii),...
        EMmodel.w(ii)*EMmodel.k, '', false)
end
xlabel('x')
ylabel('h')
title('kMeans initial estimate of map')

%Begin cycling
for ii = 1:Ncycle
    
    %draw true x meas
    x_true = rand(1,Ndraw);
    
    %true h meas
    h_true = zeros(1,Ndraw);
    for jj = 1:Ndraw
        if(x_true(jj) > 0.5)
            h_true(jj) = step_mag;
            %         pop1true = [pop1true;
            %             x_init_true(ii), h_init_true(ii)];
        else
            h_true(jj) = 0;
            %         pop2true = [pop2true;
            %             x_init_true(ii), h_init_true(ii)];
        end
    end
    
    %corrupted measurement
    meas_noise = zeros(2,Ndraw);
    for jj = 1:Ndraw
        meas_noise(:,jj) = mvnrnd([x_true(jj), h_true(jj)]', R)';
    end
    x_draw = meas_noise(1,:);
    h_draw = meas_noise(2,:);
    
    %append measurement to data
    data = [data meas_noise];
    
    %perform EM update
    disp("EM")
    tic
    EMmodel = EM(EMmodel, data, EMoptions);
    toc
    
    disp(" ")
    disp("Bayes")
    tic
    %Start Bayesian Update
    %for each measurement, find the effective number of measurements
    %belonging to each distribution
    N_eff = zeros(Ngauss,Ndraw);
    for kk = 1:Ndraw
        
        %sample all the weights
        for jj = 1:Ngauss
            N_eff(jj,kk) = gaussEval(meas_noise(:,kk),GMmodel.mu(:,jj),GMmodel.P(:,:,jj)+R);
        end
        
        %normalize weights
        N_eff(:,kk) = N_eff(:,kk)/sum(N_eff(:,kk));
        
    end
    
    %update number of measurements
    GMmodel.Nmeas = GMmodel.Nmeas + sum(N_eff,2);
    
    %now, update each element of the GM with the weighted measurment using
    %EKF
    for kk = 1:Ngauss
        
        %grab old estimates
        xbar = GMmodel.line(:,kk);
        mbar = xbar(1);
        bbar = xbar(2);
        Pbar = GMmodel.lineP(:,:,kk) + Q;
        
        %find the measurements which belong to this distribution
        targs = (N_eff(kk,:) >= 0.5);
        Ntargs = sum(targs);
        x_targs = x_draw(targs);
        h_targs = h_draw(targs);
        
        %create H, z and zbar
        H = zeros(2*Ntargs,2);
        z = zeros(2*Ntargs,1);
        zbar = zeros(2*Ntargs,1);
        Rbig = [];
        for jj = 1:Ntargs
            H(2*jj - 1,1) = -(h_targs(jj) - bbar)/mbar^2;
            H(2*jj - 1,2) = -1/mbar;
            H(2*jj, 1) = x_targs(jj);
            H(2*jj, 2) = 1;
            z(2*jj - 1) = x_targs(jj);
            z(2*jj) = h_targs(jj);
            zbar(2*jj - 1) = (h_targs(jj) - bbar)/mbar;
            zbar(2*jj) = x_targs(jj)*mbar + bbar;
            
            %update Rbig
            Rbig = blkdiag(Rbig,R);
            
        end
        
        %innovation covariance
        S = 1.5*H*Pbar*H' + Rbig;
        
        %kalman gain
        K = (Pbar*H')/S;
        
        %update line parameters
        xhat_new = xbar + K*(z - zbar);
        GMmodel.line(:,kk) = xhat_new;
        mhat = xhat_new(1);
        bhat = xhat_new(2);
        
        %covariance estimate (Josephs Form)
        I = eye(2);
        Phat_new = (I - K*H)*Pbar*(I - K*H)' + K * Rbig *K';
        GMmodel.lineP(:,:,kk) = Phat_new;
        
        %now, we must use the updated line estimate to update our estimate
        %of the mixand
        
        %for new measurements, calculate s, the distance along the line
        s_meas = zeros(Ntargs,1);
        for jj = 1:Ntargs
            %find s, this should be simplified somehow
            gamma = h_targs(jj) + x_targs(jj)/mhat;
            A = [1, -mhat;
                1, 1/mhat];
            inter = A\[bhat; gamma];
            xinter = inter(2);
            s_meas(jj) = xinter*sqrt(1+mhat^2);
            
        end
        
        %find min and max s values
        min_s = min(s_meas);
        max_s = max(s_meas);
        
        %consider updating model
        if(min_s < GMmodel.lineMin(kk))
            GMmodel.lineMin(kk) = min_s;
        else
            min_s = GMmodel.lineMin(kk);
        end
        
        if(max_s > GMmodel.lineMax(kk))
            GMmodel.lineMax(kk) = max_s;
        else
            max_s = GMmodel.lineMax(kk);
        end
        
        min_x = min_s/sqrt(1+mhat^2);
        max_x = max_s/sqrt(1+mhat^2);
        avg_x = 0.5*(min_x + max_x);
        avg_y = avg_x*mhat + bhat;
        
        %approximate the variance in the t direction
        pt_pm = avg_x/sqrt(1+mhat^2) - (bhat + mhat*avg_x - avg_y)*mhat*(1+mhat^2)^(-1.5);
        pt_pb = 1/sqrt(1+mhat^2);
        var_t = [pt_pm, pt_pb]*Phat_new*[pt_pm; pt_pb];
        
        %approximate the variance in the s direction
        k_unif = 0.75; %hueristic for scaling uniform variance
        var_s = k_unif*(1/12)*(max_s - min_s)^2;
        
        %now, rotate covariance into the cartesian frame
        P_st = [var_s, 0;
            0, var_t];
        theta = atan(mhat);
        R_st2xy = [cos(theta), -sin(theta);
            sin(theta), cos(theta)];
        P_xy = R_st2xy*P_st*R_st2xy';
        
        %update mixand
        GMmodel.mu(:,kk) = [avg_x; avg_y];
        GMmodel.P(:,:,kk) = P_xy;
        
        
    end
    
    %update total number of points
    GMmodel.Ntotal = GMmodel.Ntotal + Ndraw;
    toc
    disp(" ")
    
    
end

%sort measurements
pop1data = data(:, data(2,:) > 0.5*step_mag);
pop2data = data(:, data(2,:) <= 0.5*step_mag);

%plot elipses
for jj = 1:Ngauss
    plot_elipse(gcf, EMmodel.P(:,:,jj), EMmodel.mu(:,jj),...
        EMmodel.w(jj)*Ngauss, '', false)
    plot_elipse(gcf, GMmodel.P(:,:,jj), GMmodel.mu(:,jj),...
        GMmodel.Nmeas(jj)*Ngauss/GMmodel.Ntotal, '', false)
end
legend('Group 1','Group 2','GM 1','EM 1','GM 2','EM 2',...
    'EM 1 Final','GM 1 Final','EM 2 Final','GM 2 Final','Location','best')

figure
scatter(pop1data(1,:),pop1data(2,:))
hold on
scatter(pop2data(1,:),pop2data(2,:))
