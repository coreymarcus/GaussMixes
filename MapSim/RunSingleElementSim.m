%Mapping simulation for just one element

clear
close all
clc

%% Options

%path for sandbox
addpath("../sandbox/")
addpath("../../matlabScripts")

%truth shape
truthshape = 'Line';

%domain
x1 = -2;
x2 = 2;

%seed
rng(3);

%latex
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%estimator for line
estimator = 'KF';
% estimator = 'TLS';

%how long to stop and smell the flowers
pauselength = 1;

%% Main

%create a figure for use later
hand = figure;
hold on
plot(x1:.1:x2,TruthEval(x1:.1:x2,truthshape))

%generate some initial points
Nmeasinit = 3;
x_init = (x2 - x1)*rand(Nmeasinit,1) + x1;
y_init = TruthEval(x_init,truthshape);

%corrupt with noise
sig2 = 0.05;
x_meas = mvnrnd(x_init, sig2*eye(Nmeasinit))';
y_meas = mvnrnd(y_init, sig2*eye(Nmeasinit))';

%use TLS to initialize our estimates
[m_hat_TLS, b_hat_TLS, P_LS] = TLS(x_meas,y_meas,sig2,sig2);

%form object
obj = GaussElement(Nmeasinit);
obj.mu_mb = [m_hat_TLS; b_hat_TLS];
obj.P_mb = P_LS;

%get s1 and s2
s_meas = zeros(Nmeasinit,1);
for jj = 1:Nmeasinit
    
    %find s, this should be simplified somehow
    gamma = y_meas(jj) + x_meas(jj)/m_hat_TLS;
    A = [1, -m_hat_TLS;
        1, 1/m_hat_TLS];
    inter = A\[b_hat_TLS; gamma];
    xinter = inter(2);
    s_meas(jj) = xinter*sqrt(1+m_hat_TLS^2);
    
end
obj.s1 = min(s_meas);
obj.s2 = max(s_meas);

%create gaussians
obj = obj.Line2GaussUpdate();

%plot
% plot_handle = obj.PlotElement(hand, 1);


% Update with new measurements
Ndraw = 20;
Nupdate = 15;
for ii = 1:Nupdate
    
    %draw a new set of measurements
    x_draw = (x2 - x1)*rand(Ndraw,1) + x1;
    y_draw = TruthEval(x_draw,truthshape);
    
    %corrupt with noise
    x_meas = mvnrnd(x_draw, sig2*eye(Ndraw))';
    y_meas = mvnrnd(y_draw, sig2*eye(Ndraw))';
    
    %plot our noisy measurements
    dataplot = plot(x_meas,y_meas,'x','LineWidth',0.1,'MarkerSize',10);
    axis([-2.5 2.5 -1 6])
    
    %plotting
    plot_handle = obj.PlotElement(hand, 1);
    pause(pauselength);
    
    %delete
    delete(findall(gcf,'type','text')) %delete all exisiting text
    delete(plot_handle)
    
    %perform the update
    switch estimator
        case 'KF'
            obj = obj.UpdateLineEstimateKF(x_meas, sig2*eye(Ndraw), y_meas, sig2*eye(Ndraw));
            
        case 'TLS'
            obj = obj.UpdateLineEstimateTLS(x_meas, sig2, y_meas, sig2);
            
        otherwise
            disp('Error: invalid estimator')
    end
    
    
    %update gaussian estimate
    obj = obj.Line2GaussUpdate();
    
    %plot
    plot_handle = obj.PlotElement(hand, 1);
    axis([-2.5 2.5 -1 6])
    pause(pauselength);
    
    %delete
    if(ii ~= Nupdate)
        delete(findall(gcf,'type','text')) %delete all exisiting text
        delete(plot_handle)
        delete(dataplot)
    end
    
    
    
end

obj.mu_mb
obj.P_mb
