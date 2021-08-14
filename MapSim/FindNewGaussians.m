function [gauss_list, unused_meas] = FindNewGaussians(gauss_list, unused_meas, new_threshold, maxlength)

%Number of measurements to check
Nunuse = size(unused_meas,2);
if(Nunuse >= new_threshold)
    
    %test a variety of different clusers
    maxk = floor(Nunuse/new_threshold);
    %         testmodels = cell(maxk,1);
    %         for jj = 1:maxk
    %             testmodels{jj} = fitgmdist(unused_meas',jj,'Options',statset('MaxIter',300));
    %         end
    
    try
        model = fitgmdist(unused_meas',maxk,'Options',statset('MaxIter',300));
        modelcov = model.Sigma;
    catch
        warning('Fitgmdist failure')
        return
    end
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
            
            
        end
        
    end
    
    %remove measurements
    unused_meas = unused_meas(:,~removalidx);
    
end
end

