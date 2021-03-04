function [plot_handle, xrange, argmax_y, pdfmax] = PlotDEM(fig_handle, xhatDEM, PhatDEM, grid, x1, x2, y1, y2, dx, dy)
% plots the DEM

%get locals
xrange = x1:dx:x2;
yrange = y1:dy:y2;
Nx = length(xrange);
Ny = length(yrange);
eval_grid = zeros(Ny,Nx);

%grid
idxs = discretize(xrange,grid);

%populate
for ii = 1:Nx
    
    %gauss parameters
    mu = xhatDEM(idxs(ii));
    P = PhatDEM(idxs(ii));
    
    %find min/max idx
    [~, maxidx] = min(abs(yrange - mu - 3*sqrt(P)));
    [~, minidx] = min(abs(yrange - mu + 3*sqrt(P)));
    
    %cycle through y
    for jj = minidx:maxidx
        eval_grid(jj,ii) = gaussEval(yrange(jj),mu,P);
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

%maximum of pdf
pdfmax = max(max(eval_grid));

