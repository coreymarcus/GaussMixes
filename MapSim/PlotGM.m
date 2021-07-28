function [plot_handle, xrange, argmax_y, pdfmax] = PlotGM(fig_handle, gauss_list, x1, x2, y1, y2, d, evalmeth)
%PLOTGM plots the GM in gauss list

%get locals
Ngauss = length(gauss_list);
xrange = x1:1:x2;
yrange = y1:.05:y2;
Nx = length(xrange);
Ny = length(yrange);
eval_grid = zeros(Ny,Nx);

%cycle through gaussians
for ii = 1:Ngauss
    
    switch evalmeth
        
        case 'gauss'
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
            
        case 'mixed'
            
            %get locals
            obj = gauss_list{ii};
            Nobs = obj.Nobs;
            
            %cycle through x and y
            for jj = 1:Nx
                for kk = 1:Ny
                    
                    %eval weighted gauss at point
                    eval_pt = obj.GaussEvalST(xrange(jj),yrange(kk));
                    
                    %add
                    eval_grid(kk,jj) = eval_grid(kk,jj) + Nobs*eval_pt;
                    
                end
            end

            
            
        otherwise
            disp("Error: invalid eval method for GM!")            
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

