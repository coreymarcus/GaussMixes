function gauss_plot_handle = PlotGaussList(handle, gauss_list, meas, ...
    fitidx, truth, vehicle_state, FOV, plotsigma)

figure(handle)

Ngauss = length(gauss_list);

gauss_plot_handle = cell(Ngauss,1);

% clear exisiting data on handle
clf(handle)
hold on

% plot elements
for jj = 1:Ngauss
    gauss_plot_handle{jj} = gauss_list{jj}.PlotElement(handle, jj, plotsigma);
end

% plot measurements if needed
if(~isempty(meas))
    
    % plot unassigned measurements
    plot(meas(1,fitidx == 0), meas(2,fitidx == 0), 'x', ...
        'LineWidth',0.1,'MarkerSize',10)
    
    % plot assigned measurements
    for ii = 1:Ngauss
        plot(meas(1,fitidx == ii), meas(2,fitidx == ii), 'o', ...
            'LineWidth',0.1,'MarkerSize',10);
    end
end

% plot vehicle state if needed
if(~isempty(vehicle_state))
    
    % vehicle position
    scatter(vehicle_state(1), vehicle_state(2), 'kx')
    
    % vehicle FOV
    angle1 = vehicle_state(3) - FOV/2;
    angle2 = vehicle_state(3) + FOV/2;
    
    vec_x = [cos(angle1); cos(angle2)];
    vec_y = [sin(angle1); sin(angle2)];
    
    pos_x = [vehicle_state(1); vehicle_state(1)];
    pos_y = [vehicle_state(2); vehicle_state(2)];
    
    quiver(pos_x, pos_y, vec_x, vec_y, 2000)
    
end

% plot truth if needed
if(~isempty(truth))
    plot(truth(1,:),truth(2,:), 'k', 'LineWidth', 3)
    axis([min(truth(1,:)) max(truth(1,:)) min(truth(2,:)) max(truth(2,:))])
end



end

