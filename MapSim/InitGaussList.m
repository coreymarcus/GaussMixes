function [gauss_list, gauss_plot_handle] = InitGaussList(N_comp,meas,R)
%InitGaussList initializes a gaussian list given number of components and
%measurements

%use k-means to initialize the grouping
[~, ~, w_init, group_idx] = kMeans(N_comp, meas);

%x and y measurement
x_meas = meas(1,:);
y_meas = meas(2,:);

% form a cell array of the gaussian elements
gauss_list = cell(1,N_comp);
gauss_plot_handle = cell(1,N_comp);

%number of measurements
Nmeasinit = length(x_meas);

%cycle through the data, fitting a line to each
for ii = 1:N_comp
    
    %target idx
    targs = group_idx == ii;
    xtarg = x_meas(targs);
    ytarg = y_meas(targs);
    
    %use TLS to initialize our estimates
    [m_hat_TLS, b_hat_TLS, P_LS] = TLS(xtarg,ytarg,R(1,1),R(2,2));
    
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
    gauss_plot_handle{ii} = obj.PlotElement(gcf, ii);
    
    %list
    gauss_list{ii} = obj;
    
end

end

