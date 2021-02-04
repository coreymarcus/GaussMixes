function [plot_handle, xrange, argmax_y] = PlotGM(fig_handle, gauss_list, x1, x2, y1, y2, d)
%PLOTGM plots the GM in gauss list

%get locals
Ngauss = length(gauss_list);
xrange = x1:d:x2;
yrange = y1:d:y2;
Nx = length(xrange);
Ny = length(yrange);
eval_grid = zeros(Ny,Nx);

%cycle through gaussians
for ii = 1:Ngauss
    
    %get gauss parameters
    mu = gauss_list{ii}.mu_xy;
    P = gauss_list{ii}.P_xy;
    Nobs = gauss_list{ii}.Nobs;
    
    %cycle through x and y
    for jj = 1:Nx
        for kk = 1:Ny
            
            %eval point
            z = [xrange(jj); yrange(kk)];
            
            %eval weighted gauss at point
            eval_pt = gaussEval(z,mu,P);
            
            %add
            eval_grid(kk,jj) = eval_grid(kk,jj) + Nobs*eval_pt;
            
        end
    end
end

%find the maximum y value for each x
argmax_y = zeros(size(xrange));
for ii = 1:Nx
    [~, maxidx] = max(eval_grid(:,ii));
    argmax_y(ii) = yrange(maxidx(1));
end

%set figure
figure(fig_handle);
hold on

%plot the grid
plot_handle = mesh(xrange,yrange,eval_grid,'FaceColor','interp');
view(2)

