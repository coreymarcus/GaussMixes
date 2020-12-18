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
Ndraw = 1;

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
% EMmodel = EM(EMmodel, data, EMoptions);

%plot the initial estimate
figure
scatter(x_init_meas(group_idx == 1),h_init_meas(group_idx == 1),'xr')
hold on
scatter(x_init_meas(group_idx == 2),h_init_meas(group_idx == 2),'xb')
for ii = 1:Ngauss
    plot_elipse(gcf, P_hat_init(:,:,ii), mu_hat_init(:,ii),...
        w_hat_init(ii)*Ngauss, 'k', false)
    %     plot_elipse(gcf, EMmodel.P(:,:,ii), EMmodel.mu(:,ii),...
    %         EMmodel.w(ii)*Ngauss, 'b', false)
end
xlabel('x')
ylabel('h')
title('kMeans initial estimate of map')

%create a model for the filter
KFmodel.Nmeas = w_hat_init*Ninit;
KFmodel.Ntotal = Ninit;
KFmodel.P = P_hat_init;
KFmodel.mu = mu_hat_init;
KFmodel.line = zeros(2,Ngauss);
KFmodel.lineP = zeros(2,2,Ngauss);
KFmodel.lineMin = zeros(2,Ngauss);
KFmodel.lineMax = zeros(2,Ngauss);
KFmodel.k = Ngauss;

%use TLS to initialize the estimates of each gaussian
for ii = 1:Ngauss
    xgroup = x_init_meas(group_idx == ii);
    hgroup = h_init_meas(group_idx == ii);
    
    %call TLS
    [m, b, P] = TLS(xgroup, hgroup, P_x, P_h);
    KFmodel.line(1,ii) = m;
    KFmodel.line(2,ii) = b;
    
    %we need to project all the points to the closest intersection of the
    %line to find our "first" and "last" points
    s = zeros(length(xgroup),1);
    for jj = 1:length(xgroup)
        gamma = hgroup(jj) + xgroup(jj)/m;
        A = [1, -m;
            1, 1/m];
        inter = A\[b; gamma];
        xinter = inter(2); %yes, it is the second element
        s(jj) = xinter*sqrt(1+m^2);
    end
    
    %find min and max
    [~, minidx] = min(s);
    [~, maxidx] = max(s);
    
    KFmodel.lineMin(1,ii) = xgroup(minidx);
    KFmodel.lineMin(2,ii) = hgroup(minidx);
    KFmodel.lineMax(1,ii) = xgroup(maxidx);
    KFmodel.lineMax(2,ii) = hgroup(maxidx);
    
    %find the 
end

%Begin cycling
for ii = 1:Ncycle
    
    %draw true x meas
    x_true = rand(1,Ndraw);
    
    %true y meas
    if(x_true > 0.5)
        h_true = step_mag;
%         pop1true = [pop1true;
%             x_init_true(ii), h_init_true(ii)];
    else
        h_true = 0;
%         pop2true = [pop2true;
%             x_init_true(ii), h_init_true(ii)];
    end
    
    %corrupted measurement
    meas_noise = zeros(2,Ndraw);
    for jj = 1:Ndraw
        meas_noise(:,jj) = mvnrnd([x_true(jj), h_true(jj)]', R)';
    end
    
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
            N_eff(jj,kk) = gaussEval(meas_noise(:,kk),KFmodel.mu(:,jj),KFmodel.P(:,:,jj));
        end
        
        %normalize weights
        N_eff(:,kk) = N_eff(:,kk)/sum(N_eff(:,kk));
        
    end
    
    %now, update each element of the GM with the weighted measurment
    for kk = 1:Ngauss
        
        %update weight is a function of the effective number of
        %measurements
        w_update = zeros(2,1);
        w_update(1) = N_eff(kk);
        w_update(2) = KFmodel.Nmeas(kk);
        w_update = w_update/sum(w_update);
        
        %find new mean
        mu_update = w_update(1)*meas_noise + w_update(2)*KFmodel.mu(:,kk);
        
        %find new covariance
        P_update = w_update(1)*(R + meas_noise*meas_noise') +...
            w_update(2)*(KFmodel.P(:,:,kk) + KFmodel.mu(:,kk)*KFmodel.mu(:,kk)') - ...
            mu_update*mu_update';
        
        %update variables
        KFmodel.mu(:,kk) = mu_update;
        KFmodel.P(:,:,kk) = P_update;
        KFmodel.Nmeas(kk) = KFmodel.Nmeas(kk) + N_eff(kk);
    end
    
    %update total number of points
    KFmodel.Ntotal = KFmodel.Ntotal + Ndraw;
    toc
    disp(" ")
    
    %plot
    scatter(meas_noise(1),meas_noise(2),'rx')
    %     for jj = 1:Ngauss
    %         plot_elipse(gcf, EMmodel.P(:,:,jj), EMmodel.mu(:,jj),...
    %             EMmodel.w(jj)*Ngauss, 'b', false)
    %         plot_elipse(gcf, KFmodel.P(:,:,jj), KFmodel.mu(:,jj),...
    %             KFmodel.Nmeas(jj)*Ngauss/KFmodel.Ntotal, 'g', false)
    %     end
    
    
end

%sort measurements
pop1data = data(:, data(2,:) > 0.5*step_mag);
pop2data = data(:, data(2,:) <= 0.5*step_mag);

%plot elipses                                             
for jj = 1:Ngauss
    plot_elipse(gcf, EMmodel.P(:,:,jj), EMmodel.mu(:,jj),...
        EMmodel.w(jj)*Ngauss, 'b', false)
    plot_elipse(gcf, KFmodel.P(:,:,jj), KFmodel.mu(:,jj),...
        KFmodel.Nmeas(jj)*Ngauss/KFmodel.Ntotal, 'g', false)
end

