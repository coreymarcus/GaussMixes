function [xhat_out, Phat_out] = updateDEM(xhat, Phat, x_meas,y_meas,grid, var_y)
%updateDEM performs a basic kalman filter update on gridded points

%sort the measurments
idxs = discretize(x_meas,grid);

%number of points
N = length(xhat);

%initialize output
xhat_out = xhat;
Phat_out = Phat;

%update each one
for ii = 1:N
    
    %get points in this element
    targs = (idxs == ii);
    y_targs = y_meas(targs);
    
    %check to see if we got any
    if(isempty(y_targs))
        continue
    else
        Nupdate = length(y_targs);
    end
    
    %make sure y_targs is column
    if(~iscolumn(y_targs))
        y_targs = y_targs';
    end
    
    %get locals
    xbar = xhat(ii);
    Pbar = Phat(ii);
    
    %update
    R = var_y*eye(Nupdate);
    H = ones(Nupdate,1);
    K = Pbar*H'/(H*Pbar*H' + R);
    
    xhat_out(ii) = xbar + K*(y_targs - H*xbar);
    Phat_out(ii) = (1 - K*H)*Pbar*(1 - K*H)' + K*R*K';
    
end

end

