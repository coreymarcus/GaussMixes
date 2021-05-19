function gauss_list = UpdateGaussList(gauss_list, meas, R, fitidx, estimator)
%UpdateGaussList updates components based on chosen algo

%locals
Ngauss = length(gauss_list);
x_meas = meas(1,:);
y_meas = meas(2,:);

%perform the updates for each gaussian
for jj = 1:Ngauss
    
    %target points
    targs = fitidx == jj;
    Ntargs = sum(targs);
    
    %zero protection
    if(Ntargs <= 1)
        continue;
    end
    
    %update line estimate
    switch estimator
        case 'KF'
            gauss_list{jj} = gauss_list{jj}.UpdateLineEstimateKF(...
                x_meas(targs), R(1,1)*eye(Ntargs),...
                y_meas(targs), R(2,2)*eye(Ntargs));
            gauss_list{jj} = gauss_list{jj}.Line2GaussUpdate();
            
        case 'TLS'
            gauss_list{jj} = gauss_list{jj}.UpdateLineEstimateTLS(...
                x_meas(targs), R(1,1),...
                y_meas(targs), R(2,2));
            gauss_list{jj} = gauss_list{jj}.Line2GaussUpdate();
            
        case 'CondMerge'
            gauss_list{jj} = gauss_list{jj}.UpdateGaussDirect(...
                x_meas(targs), R(1,1), y_meas(targs), R(2,2));
            gauss_list{jj} = gauss_list{jj}.Gauss2LineUpdate();
            
        case 'NonLinLS'
            gauss_list{jj} = gauss_list{jj}.UpdateGaussNonLinLS(...
                x_meas(targs), R(1,1), y_meas(targs), R(2,2));
            gauss_list{jj} = gauss_list{jj}.Gauss2LineUpdate();
            
        case 'Direct'
            z = [x_meas(targs); y_meas(targs)];
            gauss_list{jj} = gauss_list{jj}.UpdateGaussBayes(z,R);
            gauss_list{jj} = gauss_list{jj}.Gauss2LineUpdate();
            
        otherwise
            disp('Error: invalid estimator!')
    end
    
end
end

