%Runs the mapping simulation
clear
close all
clc


%% Options

%path for sandbox
addpath("../sandbox/")
addpath("../../matlabScripts")

%truth shape
% truthshape = 'DoubleRamp'; 
truthshape = 'Parabola';
% truthshape = 'Rock';

%how do we evaluate the GM at the end
GMevalmeth = 'gauss'; %uses the standard gaussian
% GMevalmeth = 'mixed'; %uses a mixed representation based on s and t

% how do we choose our slope estimate?
slopemethod = "ML"; %maximum liklihood
% slopemethod = "MMSE";

%how long to stop and smell the flowers
pauselength = 0;

%maximum lenght of an element
maxlength = 1;

%threshold for merging two gaussians
mergethresh = 0.25;

%estimator for line
% estimator = 'KF';
% estimator = 'TLS';
% estimator = 'CondMerge';
% estimator = 'NonLinLS'; %nonlin LS
estimator = 'Direct'; %direct estimation of the gaussian


%domain
% x1 = -2;
% x2 = 2;
x1 = -10;
x2 = 10;
% Ndem = floor((x2 - x1)/maxlength)+1; %number of DEM bins
Ndem = 10; %number of DEM bins

%seed
rng(3);

%measurement noise covariance
sig2 = 0.25;
% sig2 = 0.01;
% sig2 = 0;

%latex
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Main

%generate some initial points
Nmeasinit = 50;
% x_init = 0.3*(x2 - x1)*rand(Nmeasinit,1) + randi([0,1],Nmeasinit,1)*x1;
x_init = (x2 - x1)*rand(Nmeasinit,1) + x1;
y_init = TruthEval(x_init,truthshape);

%corrupt with noise
x_meas = mvnrnd(x_init, sig2*eye(Nmeasinit))';
y_meas = mvnrnd(y_init, sig2*eye(Nmeasinit))';

%use k-means to initialize the grouping
Ngauss = 2;
[~, ~, w_init, group_idx] = kMeans(Ngauss, [x_meas, y_meas]');

%create a figure for use later
hand = figure;
plot(x1:.1:x2,TruthEval(x1:.1:x2,truthshape))

%initialize the DEM estimate
xhatDEM = zeros(Ndem,1);
PhatDEM = ones(Ndem,1);
grid = linspace(x1-1,x2+1,Ndem);
[xhatDEM, PhatDEM] = updateDEM(xhatDEM, PhatDEM, x_meas, y_meas, grid, sig2);

%pool of un-used measurements
unused_meas = [];

% form a cell array of the gaussian elements
gauss_list = cell(1,Ngauss);
gauss_plot_handle = cell(1,Ngauss);
for ii = 1:Ngauss
    
    %target idx
    targs = group_idx == ii;
    xtarg = x_meas(targs);
    ytarg = y_meas(targs);
    
    %use TLS to initialize our estimates
    [m_hat_TLS, b_hat_TLS, P_LS] = TLS(xtarg,ytarg,sig2,sig2);
    
    %form object
    obj = GaussElement(w_init(ii)*Nmeasinit);
    obj.mu_mb = [m_hat_TLS; b_hat_TLS];
    obj.P_mb = P_LS;
    
    %get s1 and s2
    s_meas = zeros(length(xtarg),1);
    for jj = 1:length(xtarg)
        
        %find s, this should be simplified somehow
        gamma = ytarg(jj) + xtarg(jj)/m_hat_TLS;
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
    
    %initialize parameters for bayesian inference
    obj.n_dof = 4;
    obj.Psi = (obj.n_dof - 2  - 1)*obj.P_xy;
    
    %plot
    gauss_plot_handle{ii} = obj.PlotElement(hand, ii);
    
    %list
    gauss_list{ii} = obj;
    
end

% Update with new measurements
Ndraw = 50;
Nupdate = 5;
for ii = 1:Nupdate
    
    %draw a new set of measurements
    %x_draw = 0.3*(x2 - x1)*rand(Ndraw,1) + randi([0,1],Ndraw,1)*x1;
    x_draw = (x2 - x1)*rand(Ndraw,1) + x1;
    y_draw = TruthEval(x_draw,truthshape);
    
    %corrupt with noise
    x_meas = mvnrnd(x_draw, sig2*eye(Ndraw))';
    y_meas = mvnrnd(y_draw, sig2*eye(Ndraw))';
    
    %update the DEM
    [xhatDEM, PhatDEM] = updateDEM(xhatDEM, PhatDEM, x_meas, y_meas, grid, sig2);
    
    
    %     %evaluate the liklihood that each measurement came from each gaussian
    %     wdraw = zeros(Ndraw,Ngauss);
    %     fitidx = zeros(Ndraw,1);
    %     for jj = 1:Ndraw
    %         for kk = 1:Ngauss
    %             wdraw(jj,kk) = gauss_list{kk}.Nobs*gauss_list{kk}.GaussEvalST(x_meas(jj),y_meas(jj),sig2*eye(2));
    %         end
    %
    %         %determine which one is the best fit
    %         [maxval, maxidx] = max(wdraw(jj,:));
    %         if(maxval > 0)
    %             fitidx(jj) = maxidx(1);
    %         end
    %
    %     end
    
    %evaluate the sigma distance of each measurement to each gaussian
    sigdist = zeros(Ndraw,Ngauss);
    fitidx = zeros(Ndraw,1);
    for jj = 1:Ndraw
        for kk = 1:Ngauss
            sigdist(jj,kk) = gauss_list{kk}.GaussSig(x_meas(jj),y_meas(jj),sig2*eye(2));
        end
        
        %determine which one is the best fit
        [minval, minidx] = min(sigdist(jj,:));
        if(minval < 3)
            fitidx(jj) = minidx(1);
        else %% add to list of unused measurements
            unused_meas = [unused_meas, [x_meas(jj); y_meas(jj)]];            
        end
    end
    
    %plot our noisy measurements
    dataplot = plot(x_meas(fitidx == 0),y_meas(fitidx == 0),'x',...
        x_meas(fitidx == 1),y_meas(fitidx == 1),'o',...
        x_meas(fitidx == 2),y_meas(fitidx == 2),'o',...
        x_meas(fitidx == 3),y_meas(fitidx == 3),'o',...
        x_meas(fitidx == 4),y_meas(fitidx == 4),'o',...
        x_meas(fitidx == 5),y_meas(fitidx == 5),'o',...
        x_meas(fitidx == 6),y_meas(fitidx == 6),'o',...
        'LineWidth',0.1,'MarkerSize',10);
    axis([1.2*x1 1.2*x2 1.2*x1 1.2*x2])
    pause(pauselength);
    
    %plotting
    for jj = 1:Ngauss
        delete(gauss_plot_handle{jj})
        
        gauss_plot_handle{jj} = gauss_list{jj}.PlotElement(hand, jj);
    end
    
    %perform the updates for each gaussian
    gauss_list_new = cell(0);
    gauss_plot_handle_new = cell(0);
    for jj = 1:Ngauss
        
        %target points
        targs = fitidx == jj;
        Ntargs = sum(targs);
        
        %zero protection
        if(Ntargs <= 1)
            continue;
        end
%         
%         if(ii == 7 && jj == 5)
%             disp('pause here')
%         end
        
        %update line estimate
        switch estimator
            case 'KF'
                gauss_list{jj} = gauss_list{jj}.UpdateLineEstimateKF(...
                    x_meas(targs), sig2*eye(Ntargs),...
                    y_meas(targs), sig2*eye(Ntargs));
                gauss_list{jj} = gauss_list{jj}.Line2GaussUpdate();
                
            case 'TLS'
                gauss_list{jj} = gauss_list{jj}.UpdateLineEstimateTLS(...
                    x_meas(targs), sig2,...
                    y_meas(targs), sig2);
                gauss_list{jj} = gauss_list{jj}.Line2GaussUpdate();
                
            case 'CondMerge'
                gauss_list{jj} = gauss_list{jj}.UpdateGaussDirect(...
                    x_meas(targs), sig2, y_meas(targs), sig2);
                gauss_list{jj} = gauss_list{jj}.Gauss2LineUpdate();
                
            case 'NonLinLS'
                gauss_list{jj} = gauss_list{jj}.UpdateGaussNonLinLS(...
                    x_meas(targs), sig2, y_meas(targs), sig2);
                gauss_list{jj} = gauss_list{jj}.Gauss2LineUpdate();
                
            case 'Direct'
                z = [x_meas(targs)'; y_meas(targs)'];
                R = sig2*eye(2);
                gauss_list{jj} = gauss_list{jj}.UpdateGaussBayes(z,R);
                gauss_list{jj} = gauss_list{jj}.Gauss2LineUpdate();
                
            otherwise
                disp('Error: invalid estimator!')
        end
        
        %check to see if we should split this gaussian
        if((gauss_list{jj}.s2 - gauss_list{jj}.s1) > maxlength)
            [gauss_list{jj}, gauss_list_new{end+1}] = gauss_list{jj}.SplitElement();
            gauss_plot_handle_new{end+1} = gauss_plot_handle{jj};
        end
        
    end
    
    %update list of gaussian
    gauss_list = [gauss_list, gauss_list_new];
    gauss_plot_handle = [gauss_plot_handle, gauss_plot_handle_new];
    
    %update number of gaussians
    Ngauss = length(gauss_list);
    
    %plotting
    delete(findall(gcf,'type','text')) %delete all exisiting text
    for jj = 1:Ngauss
        delete(gauss_plot_handle{jj})
        
        gauss_plot_handle{jj} = gauss_list{jj}.PlotElement(hand, jj);
    end
    
    %check to see if we can discover new gaussians from any of the unused
    %measurements
    
    Nunuse = size(unused_meas,2);
    if(Nunuse >= 40)
        
        %test a variety of different clusers
        maxk = floor(Nunuse/40);
        %         testmodels = cell(maxk,1);
        %         for jj = 1:maxk
        %             testmodels{jj} = fitgmdist(unused_meas',jj,'Options',statset('MaxIter',300));
        %         end
        
        model = fitgmdist(unused_meas',maxk,'Options',statset('MaxIter',300));
        modelcov = model.Sigma;
        
        %initialize idxs for removal
        removalidx = zeros(Nunuse,1);
        
        %look at each element, check to see if it is small enough to be
        %accepted
        for jj = 1:maxk
            
            %eigen decompose
            D = eig(modelcov(:,:,jj));

            if(max(D) < 0.5*maxlength)
                %this should be a new element!
                
                %cluster the data
                clusteridx = cluster(model,unused_meas');
                
                %get the target data
                targidx = (clusteridx == jj);
                removalidx = removalidx + targidx;
                
                %form object
                obj = GaussElement(sum(targidx));
                obj.P_xy = modelcov(:,:,jj);
                obj.mu_xy = model.mu(jj,:)';
                obj.Nobs = sum(targidx);
                obj = obj.Gauss2LineUpdate();
                
                %initialize parameters for bayesian inference
                obj.n_dof = 4;
                obj.Psi = (obj.n_dof - 2 - 1)*obj.P_xy;
                
                % add
                gauss_list{end+1} = obj;
                gauss_plot_handle{end+1} = gauss_list{end}.PlotElement(hand,0);
                
                
            end
            
        end
        
        %remove measurements
        unused_meas = unused_meas(:,~removalidx); 
            
    end
    
    %update number of gaussians
    Ngauss = length(gauss_list);
       
    %pause for inspection
    axis([1.2*x1 1.2*x2 1.2*x1 1.2*x2])
    pause(pauselength);
    delete(dataplot);
    
    %check the sigma length between gaussians for merging
    sigdist = zeros(Ngauss);
    for jj = 1:(Ngauss-1)
        for kk = (jj+1):Ngauss
            
            %find sigma distance between the two gaussians
            mu_xy = gauss_list{jj}.mu_xy;
            P_xy = gauss_list{jj}.P_xy;
            sigdist(jj,kk) = gauss_list{kk}.GaussSig(mu_xy(1),mu_xy(2),P_xy);
        end
    end
    
    %see if any of them should be merged
    mergedidxs = [];
    for jj = 1:(Ngauss - 1)
        
        %extract row
        rowvals = sigdist(jj,:);
        rowvals(rowvals <= 0) = inf;
        
        %find minimum
        [minval, minidx] = min(rowvals);
        
        %if we are less than some threshold, merge
        if(minval <= mergethresh)
            
            %we will get rid of the jj'th element
            mergedidxs = [mergedidxs jj];
            
            %extract objects
            obj1 = gauss_list{jj};
            obj2 = gauss_list{minidx};
            
            %locals
            N1 = obj1.Nobs;
            N2 = obj2.Nobs;
            w1 = N1/(N1 + N2);
            w2 = N2/(N1 + N2);
            
            %new object
            obj3 = GaussElement(N1 + N2);
            obj3.mu_xy = w1*obj1.mu_xy + w2*obj2.mu_xy;
            obj3.P_xy = w1*(obj1.P_xy + obj1.mu_xy*obj1.mu_xy')...
                + w2*(obj2.P_xy + obj2.mu_xy*obj2.mu_xy')...
                -obj3.mu_xy*obj3.mu_xy';
            obj3 = obj3.Gauss2LineUpdate();
            
            %initialize parameters for bayesian inference
            obj3.n_dof = 4;
            obj3.Psi = (obj3.n_dof - 2  - 1)*obj3.P_xy;
            
            %assign
            gauss_list{minidx} = obj3;
            
        end
    end
    
    %get rid of merged objects
    toremove = zeros(Ngauss,1);
    toremove(mergedidxs) = 1;
    gauss_list = gauss_list(~toremove);
    
    %get rid of all the plots
    delete(findall(gcf,'type','text')) %delete all exisiting text
    for jj = 1:Ngauss
        delete(gauss_plot_handle{jj})
    end
    
    %remove merged plots
    gauss_plot_handle = gauss_plot_handle(~toremove);
    
    %update number of gaussians
    Ngauss = length(gauss_list);
    
    %plot
    for jj = 1:Ngauss
        gauss_plot_handle{jj} = gauss_list{jj}.PlotElement(hand, jj);
    end
    pause(pauselength)

end


%% Plotting

%plot the final results
pdf_plot = figure;
[~, xmax, yargmax, maxGM] = PlotGM(pdf_plot,gauss_list,x1,x2,x1,x2,.05,GMevalmeth);
plot3(x1:.1:x2,TruthEval(x1:.1:x2,truthshape),maxGM*ones(length(x1:.1:x2),1),'k','LineWidth',2)
plot3(xmax, yargmax,maxGM*ones(length(xmax),1),'r','LineWidth',2)
xlabel('x')
ylabel('y')
legend('PDF','True Surface','Estimated Surface')
title('PDF From GM')

DEM_plot = figure;
[~, xmaxDEM, yargmaxDEM, maxDEM] = PlotDEM(DEM_plot, xhatDEM, PhatDEM, grid, x1,x2,x1,x2,.1,.001);
plot3(x1:.1:x2,TruthEval(x1:.1:x2,truthshape),maxDEM*ones(length(x1:.1:x2),1),'k','LineWidth',2)
plot3(xmaxDEM, yargmaxDEM, maxDEM*ones(length(xmaxDEM),1),'r','LineWidth',2)
title('PDF from DEM')


%get the slope estimate at each point
xsample = x1:.05:x2;
Nsample = length(xsample);
slopeGM = zeros(Nsample,1);
slopeDEM = zeros(Nsample,1);
for ii = 1:Nsample
    slopeGM(ii) = GetGMSlopeEstimate(gauss_list,xsample(ii),slopemethod);
    slopeDEM(ii) = GetDEMSlopeEstimate(xhatDEM, grid, xsample(ii));
end
[~, slopetruth] = TruthEval(xsample,truthshape);

%plot
GMslope_plot = figure;
plot(xsample,slopeGM,'LineWidth',2);
hold on
plot(xsample,slopetruth,'LineWidth',2);
title({'GM Slope Estimate',slopemethod})

%plot
DEMslope_plot = figure;
plot(xsample,slopeDEM,'LineWidth',2);
hold on
plot(xsample,slopetruth,'LineWidth',2);
title('DEM slope estimate')