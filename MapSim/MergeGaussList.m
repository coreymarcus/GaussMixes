function gauss_list = MergeGaussList(gauss_list, mergethreshold)
%UpdateGaussList updates components based on chosen algo

%locals
Ngauss = length(gauss_list);

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
    if(minval <= mergethreshold)
        
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


end

