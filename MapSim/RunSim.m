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


%% Main

%generate some initial points
x1 = -10;
x2 = 10;
Nmeasinit = 50;
% x_init = 0.3*(x2 - x1)*rand(Nmeasinit,1) + randi([0,1],Nmeasinit,1)*x1;
x_init = (x2 - x1)*rand(Nmeasinit,1) + x1;
y_init = TruthEval(x_init,truthshape);

%corrupt with noise
sig2 = 0.1;
x_meas = mvnrnd(x_init, sig2*eye(Nmeasinit))';
y_meas = mvnrnd(y_init, sig2*eye(Nmeasinit))';

%use k-means to initialize the grouping
Ngauss = 5;
[~, ~, w_init, group_idx] = kMeans(Ngauss, [x_meas, y_meas]');

%create a figure for use later
hand = figure;
plot(x1:.1:x2,TruthEval(x1:.1:x2,truthshape))

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
    
    %plot
    gauss_plot_handle{ii} = obj.PlotElement(hand);
    
    %list
    gauss_list{ii} = obj;
    
end

% Update with new measurements
Ndraw = 50;
Nupdate = 15;
for ii = 1:Nupdate
    
    %draw a new set of measurements
    %x_draw = 0.3*(x2 - x1)*rand(Ndraw,1) + randi([0,1],Ndraw,1)*x1;
    x_draw = (x2 - x1)*rand(Ndraw,1) + x1;
    y_draw = TruthEval(x_draw,truthshape);
    
    %corrupt with noise
    x_meas = mvnrnd(x_draw, sig2*eye(Ndraw))';
    y_meas = mvnrnd(y_draw, sig2*eye(Ndraw))';
    
    
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
        if(minidx < 2)
            fitidx(jj) = minidx(1);
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
    pause(1);
    
    %plotting
    for jj = 1:Ngauss
        delete(gauss_plot_handle{jj})
        
        gauss_plot_handle{jj} = gauss_list{jj}.PlotElement(hand);
    end
    
    %perform the updates for each gaussian
    gauss_list_new = cell(0);
    gauss_plot_handle_new = cell(0);
    for jj = 1:Ngauss
        
        %target points
        targs = fitidx == jj;
        Ntargs = sum(targs);
        
        %zero protection
        if(Ntargs == 0)
            continue;
        end
        
        %update line estimate
        gauss_list{jj} = gauss_list{jj}.UpdateLineEstimateKF(...
            x_meas(targs), sig2*eye(Ntargs),...
            y_meas(targs), sig2*eye(Ntargs));
        
        %update gaussian estimate
        gauss_list{jj} = gauss_list{jj}.Line2GaussUpdate();
        
        %check to see if we should split this gaussian
        if((gauss_list{jj}.s2 - gauss_list{jj}.s1) > 1)
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
    for jj = 1:Ngauss
        delete(gauss_plot_handle{jj})
        
        gauss_plot_handle{jj} = gauss_list{jj}.PlotElement(hand);
    end
    
    %pause for inspection
    axis([1.2*x1 1.2*x2 1.2*x1 1.2*x2])
    pause(1);
    delete(dataplot);
    
    
    
end

%plot the final results
pdf_plot = figure;
[~, xmax, yargmax] = PlotGM(pdf_plot,gauss_list,x1,x2,x1,x2,.05);
plot3(x1:.1:x2,TruthEval(x1:.1:x2,truthshape),100*ones(length(x1:.1:x2),1),'k','LineWidth',2)
plot3(xmax, yargmax,100*ones(length(xmax),1),'r','LineWidth',2)