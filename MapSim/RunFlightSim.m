clear
close all
clc


%% Options

%path for sandbox
addpath("../sandbox/")
addpath("../../matlabScripts")

%how do we evaluate the GM at the end?
GMevalmeth = 'gauss'; %uses the standard gaussian
% GMevalmeth = 'mixed'; %uses a mixed representation based on s and t

%create plots during execution?
statusplots = true;

%time parameters
t_sim = 59;
t_step = 1;
pause_length = 2; % time to pause and look at plots for

%trajectory
traj = 'Ramp'; %straight ramp descent

% true terrain
terrain = 'GentleParabola';

% measurement model
meas_model = 'Simple';

%initialization settings
N_meas_init = 1000;
R_init = eye(2);
N_gauss_init = 3;

%cycle settings
N_meas_cycle = 50;
FOV = pi/6;
R_cycle = eye(2);

%gaussian parameters
maxlength = 10;
mergethreshold = 0.08;
bindist = 1.5;
plotsigma = 2*sqrt(mergethreshold);

%estimator to use
estimator = 'KF';
% estimator = 'TLS';
% estimator = 'CondMerge';
% estimator = 'NonLinLS'; %nonlin LS
% estimator = 'Direct'; %direct estimation of the gaussian

%% Main

%figure for plotting
main_fig = figure;

%draw initial measurements
x1 = -400;
x2 = 400;
x_init = (x2 - x1)*rand(N_meas_init,1) + x1;
y_init = TruthEval(x_init,terrain);

%corrupt with noise
x_meas = mvnrnd(x_init, R_init(1,1)*eye(N_meas_init))';
y_meas = mvnrnd(y_init, R_init(2,2)*eye(N_meas_init))';

%initialize list
[gauss_list, gauss_plot_handle] = InitGaussList(N_gauss_init, [x_meas, y_meas]', R_init);
for jj = 1:N_gauss_init
    delete(gauss_plot_handle{jj})
end

%sample truth
truth_samp = [x1:.1:x2;
    TruthEval(x1:.1:x2,terrain)];

%time vector
t = 0:t_step:t_sim;
N_t = length(t);

%initialize storage
x = zeros(3,N_t);

%initialize measurement storage vector
meas = zeros(2,N_t*N_meas_cycle);

%theta values for measurement FOV
theta = linspace(-FOV/2, FOV/2, N_meas_cycle);

%initialize unused measurements list
unused_meas = [];

%loop
for ii = 1:N_t
    
    %get true state
    x_iter = GetVehiclePose(t(ii),traj);
    
    %measurement noise
    noise_iter = mvnrnd([0; 0],R_cycle,N_meas_cycle)';
    
    %get measurment
    meas_iter = zeros(2,N_meas_cycle);
    for jj = 1:N_meas_cycle
        meas_iter(:,jj) = GetTerrainMeasurement(x_iter,theta(jj),zeros(2),meas_model,terrain) ...
            + noise_iter(:,jj);
    end
    
    %bin the measurements
    fitidx = BinMeasurements(meas_iter,gauss_list,R_cycle, bindist);
    
    %build list of unused measurements
    unused_meas = [unused_meas, meas_iter(:,fitidx == 0)];
    
    %look for potential new gaussians
    [gauss_list, unused_meas] = FindNewGaussians(gauss_list, unused_meas, ...
        40, maxlength);
    
    %update the list
    gauss_list = UpdateGaussList(gauss_list, meas_iter, R_cycle, fitidx, ...
        estimator, maxlength);
    
    % plot
%     if(statusplots)
%         gauss_plot_handle = PlotGaussList(main_fig, gauss_list, meas_iter, fitidx, truth_samp, ...
%             x_iter, FOV, plotsigma);
%         pause(pause_length);
%     end
    
    %consider merging
    %     gauss_list = MergeGaussList(gauss_list, mergethreshold);
    
    % plot
    if(statusplots)
        gauss_plot_handle = PlotGaussList(main_fig, gauss_list, meas_iter, fitidx, truth_samp, ...
            x_iter, FOV, plotsigma);
        pause(pause_length)
    end
    
    %store
    x(:,ii) = x_iter;
    meas(:,(ii-1)*N_meas_cycle+1:ii*N_meas_cycle) = meas_iter;
    
    disp(ii/N_t)
    
end

%% Plotting

figure
plot(x(1,:),x(2,:),'*')
axis equal
hold on
plot(truth_samp(1,:),truth_samp(2,:))
xlabel('x')
ylabel('y')
legend('Vehicle Trajectory','Terrain')

%plot the final results
pdf_plot = figure;
[~, xmax, yargmax, maxGM] = PlotGM(pdf_plot,gauss_list,x1,x2,-20,40,.1,GMevalmeth);
plot3(truth_samp(1,:),truth_samp(2,:),maxGM*ones(size(truth_samp,2),1),'k','LineWidth',2)
plot3(xmax, yargmax,maxGM*ones(length(xmax),1),'r','LineWidth',2)
xlabel('x')
ylabel('y')
legend('PDF','True Surface','Estimated Surface')
title('PDF From GM')

%plot
final_plot = figure;
for jj = 1:length(gauss_list)
    gauss_plot_handle{jj} = gauss_list{jj}.PlotElement(final_plot, jj);
end
plot(x1:1:x2,TruthEval(x1:1:x2,terrain))

figure
subplot(2,1,1)
plot(t,x(1,:))

subplot(2,1,2)
plot(t,x(2,:))



