classdef PlanarCell
    %Planar Variable resolution 3D DEM cell tree structure
    
    properties
        isDivided_ = false; % Is this cell divided?
        isValid_ = false; % Does this cell have a valid estiamte?
        northwest_; % Child cell if divided
        northeast_; % Child cell if divided
        southwest_;
        southeast_;
        center_; % [2x1] cell center in global frame
        halfWidth_;
        halfHeight_;
        nhat_; % Surface normal estimate
        dhat_; % Orthogonal distance from origin estimate
        Phat_; % [3x3] error covariance estimate
        points_; %[3xN] array of measurements in the global frame
        Rpoints_; %[3x3xN] array of measurement covariances in global frame
        nPointsMax_; % Number of observations after which division occurs
        depth_; % Depth of this cell ( = 0 for top level, = positive for subsequent levels)
        maxDepth_; % Maximum depth for division
        numOutliers_ = 0; % Number of outliers detected
        numReject_ = 0; % Number of measurements rejected
        rejectThresh_ = 500; % Cell reinitialized if more than this number of measurements rejected
        fitScore_ = 0; % Sum of outlier residuals
        avgFitScoreThresh_ = 1.5; % Average fit score required for division
    end
    
    methods
        
        function obj = AddPoints(obj, rmat, Rmat, bhatmat, m0mat)
            
            %If this cell is divided, add points to children
            if(obj.isDivided_)
                % Divide the points
                xbins = discretize(rmat(1,:),[obj.center_(1) - obj.halfWidth_,...
                    obj.center_(1),...
                    obj.center_(1) + obj.halfWidth_]);
                ybins = discretize(rmat(2,:),[obj.center_(2) - obj.halfHeight_,...
                    obj.center_(2),...
                    obj.center_(2) + obj.halfHeight_]);
                nwpoints = xbins == 1 & ybins == 2;
                nepoints = xbins == 2 & ybins == 2;
                swpoints = xbins == 1 & ybins == 1;
                sepoints = xbins == 2 & ybins == 1;
                
                npts = size(bhatmat,2);
                nptsactual = sum(nwpoints) + sum(nepoints) + sum(swpoints) + sum(sepoints);
                if(npts ~= nptsactual)
                    warning("Losing points!")
                end
                
                % Add points to each child
                obj.northwest_ = obj.northwest_.AddPoints(rmat(:,nwpoints),...
                    Rmat(:,:,nwpoints),...
                    bhatmat(:,nwpoints),...
                    m0mat(:,nwpoints));
                obj.northeast_ = obj.northeast_.AddPoints(rmat(:,nepoints),...
                    Rmat(:,:,nepoints),...
                    bhatmat(:,nepoints),...
                    m0mat(:,nepoints));
                obj.southwest_ = obj.southwest_.AddPoints(rmat(:,swpoints),...
                    Rmat(:,:,swpoints),...
                    bhatmat(:,swpoints),...
                    m0mat(:,swpoints));
                obj.southeast_ = obj.southeast_.AddPoints(rmat(:,sepoints),...
                    Rmat(:,:,sepoints),...
                    bhatmat(:,sepoints),...
                    m0mat(:,sepoints));
                return
            end
            
            
            % Locals
            nmeas = size(rmat,2);
            
            % Potential total number of points in cell
            npoints = size(obj.points_,2) + size(rmat,2);
            
            % Check for outliers
            if(obj.isValid_)
                 obj = obj.CountOutliers(rmat, Rmat, bhatmat, m0mat);
            end
            
            dividelogic = npoints > 30 & ...
                obj.fitScore_/npoints > obj.avgFitScoreThresh_ & ...
                obj.depth_ < obj.maxDepth_;
            
            % If the number of points is greater than the max, divide
            if(dividelogic)
            %if(npoints > obj.nPointsMax_ && obj.depth_ < obj.maxDepth_)
                
                % add points to points
                if(size(obj.points_,2) > 0)
                    obj.points_(:,end+1:end+nmeas) = rmat;
                    obj.Rpoints_(:,:,end+1:end+nmeas) = Rmat;
                else
                    obj.points_ = rmat;
                    obj.Rpoints_ = Rmat;
                end
                
                obj = obj.Divide();
                
                return
            end
            
            % Either initialize or update the cell
            if(~obj.isValid_)
                
                % add points to points
                if(size(obj.points_,2) > 0)
                    obj.points_(:,end+1:end+nmeas) = rmat;
                    obj.Rpoints_(:,:,end+1:end+nmeas) = Rmat;
                else
                    obj.points_ = rmat;
                    obj.Rpoints_ = Rmat;
                end
                
                % initialize
                obj = obj.InitFromML();
            else
                % Cycle, performing updates
                for ii = 1:nmeas
                    obj = obj.UpdateEKF(rmat(:,ii),...
                        Rmat(:,:,ii),...
                        bhatmat(:,ii),...
                        m0mat(:,ii));
                end
                
                % Check the total number of rejected measurements
                if obj.numReject_ > obj.rejectThresh_
                    % Reinitialize
                    obj = obj.InitFromML();
                    warning('Cell Reinitialized')                    
                end
            end
        end
        
        function obj = UpdateEKF(obj, r, R, bhat, m0)
            
            % Locals
            nhat = obj.nhat_;
            dhat = obj.dhat_;
            Pbar = obj.Phat_;
            
            % Transform measurement to cell frame
            rcell = r;
            rcell(1:2) = rcell(1:2) - obj.center_;
            m0cell = m0;
            m0cell(1:2) = m0cell(1:2) - obj.center_;
            
            % Predicted measurement
            t_exp = (dhat - nhat'*m0cell)/(nhat'*bhat);
            r_exp = m0cell + t_exp*bhat;
            
            % Skip if observation angle is too large or expected
            % range negative
            theta = acos(nhat'*bhat)*180/pi;
            if(theta < 110 || t_exp <= 0)
                warning('Measurement rejected')
                obj.numReject_ = obj.numReject_ + 1;
                obj.points_(:,end+1) = r;
                obj.Rpoints_(:,:,end+1) = R;
                return;
            end
            
            % Find the measurement jacobian
            H = obj.JacobianEval(bhat(1),bhat(2),bhat(3),dhat, m0cell(1),m0cell(2),m0cell(3),nhat(1),nhat(2),nhat(3));
            
            % Find the kalman gain
            K = Pbar*H'/(H*Pbar*H' + R);
            
            % Find the update
            update = K*(rcell - r_exp);
            
            % Perform the update
            nhatplus = UnitVectorAdd(nhat,update(1:2));
            dhatplus = dhat + update(3);
            Phat = (eye(3) - K*H)*Pbar*(eye(3) - K*H)' + K*R*K';
            
            % Assign the update
            obj.nhat_ = nhatplus;
            obj.dhat_ = dhatplus;
            obj.Phat_ = Phat;
            obj.points_(:,end+1) = r;
            obj.Rpoints_(:,:,end+1) = R;
        end
        
        function J = JacobianEval(~,bhat_1,bhat_2,bhat_3,d,m0_1,m0_2,m0_3,nhat_1,nhat_2,nhat_3)
            %JACOBIANEVAL
            %    J = JACOBIANEVAL(BHAT_1,BHAT_2,BHAT_3,D,M0_1,M0_2,M0_3,NHAT_1,NHAT_2,NHAT_3)
            
            %    This function was generated by the Symbolic Math Toolbox version 8.1.
            %    11-Aug-2021 19:28:24
            
            t2 = nhat_3+1.0;
            t3 = 1.0./t2;
            t4 = nhat_1.^2;
            t5 = t3.*t4;
            t6 = t5-1.0;
            t7 = bhat_1.*nhat_1;
            t8 = bhat_2.*nhat_2;
            t9 = bhat_3.*nhat_3;
            t10 = t7+t8+t9;
            t11 = 1.0./t10;
            t12 = nhat_2.^2;
            t13 = t3.*t12;
            t14 = t13-1.0;
            t15 = 1.0./t10.^2;
            t16 = m0_1.*nhat_1;
            t17 = m0_2.*nhat_2;
            t18 = m0_3.*nhat_3;
            t19 = -d+t16+t17+t18;
            t20 = m0_3.*nhat_1;
            t21 = m0_1.*t6;
            t22 = m0_2.*nhat_1.*nhat_2.*t3;
            t23 = t20+t21+t22;
            t24 = bhat_3.*nhat_1;
            t25 = bhat_1.*t6;
            t26 = bhat_2.*nhat_1.*nhat_2.*t3;
            t27 = t24+t25+t26;
            t28 = m0_3.*nhat_2;
            t29 = m0_2.*t14;
            t30 = m0_1.*nhat_1.*nhat_2.*t3;
            t31 = t28+t29+t30;
            t32 = bhat_3.*nhat_2;
            t33 = bhat_2.*t14;
            t34 = bhat_1.*nhat_1.*nhat_2.*t3;
            t35 = t32+t33+t34;
            J = reshape([bhat_1.*t11.*t23-bhat_1.*t15.*t19.*t27,bhat_2.*t11.*t23-bhat_2.*t15.*t19.*t27,bhat_3.*t11.*t23-bhat_3.*t15.*t19.*t27,bhat_1.*t11.*t31-bhat_1.*t15.*t19.*t35,bhat_2.*t11.*t31-bhat_2.*t15.*t19.*t35,bhat_3.*t11.*t31-bhat_3.*t15.*t19.*t35,bhat_1.*t11,bhat_2.*t11,bhat_3.*t11],[3,3]);
            
        end
        
        function obj = Divide(obj)
            
            % initialize children
            northwest = PlanarCell();
            northeast = PlanarCell();
            southwest = PlanarCell();
            southeast = PlanarCell();
            
            % Assign children bounds
            northwest.center_ = obj.center_ + [-obj.halfWidth_;
                obj.halfHeight_]/2;
            northeast.center_ = obj.center_ + [obj.halfWidth_;
                obj.halfHeight_]/2;
            southwest.center_ = obj.center_ + [-obj.halfWidth_;
                -obj.halfHeight_]/2;
            southeast.center_ = obj.center_ + [obj.halfWidth_;
                -obj.halfHeight_]/2;
            northwest.halfWidth_ = obj.halfWidth_/2;
            northeast.halfWidth_ = obj.halfWidth_/2;
            southwest.halfWidth_ = obj.halfWidth_/2;
            southeast.halfWidth_ = obj.halfWidth_/2;
            northwest.halfHeight_ = obj.halfHeight_/2;
            northeast.halfHeight_ = obj.halfHeight_/2;
            southwest.halfHeight_ = obj.halfHeight_/2;
            southeast.halfHeight_ = obj.halfHeight_/2;
            
            % assign children depth
            northwest.depth_ = obj.depth_ + 1;
            northeast.depth_ = obj.depth_ + 1;
            southwest.depth_ = obj.depth_ + 1;
            southeast.depth_ = obj.depth_ + 1;
            
            % assign children other parameters
            northwest.nPointsMax_ = obj.nPointsMax_;
            northeast.nPointsMax_ = obj.nPointsMax_;
            southwest.nPointsMax_ = obj.nPointsMax_;
            southeast.nPointsMax_ = obj.nPointsMax_;
            northwest.maxDepth_ = obj.maxDepth_;
            northeast.maxDepth_ = obj.maxDepth_;
            southwest.maxDepth_ = obj.maxDepth_;
            southeast.maxDepth_ = obj.maxDepth_;
            
            % Divide the points
            xbins = discretize(obj.points_(1,:),[obj.center_(1) - obj.halfWidth_,...
                obj.center_(1),...
                obj.center_(1) + obj.halfWidth_]);
            ybins = discretize(obj.points_(2,:),[obj.center_(2) - obj.halfHeight_,...
                obj.center_(2),...
                obj.center_(2) + obj.halfHeight_]);
            nwpoints = xbins == 1 & ybins == 2;
            nepoints = xbins == 2 & ybins == 2;
            swpoints = xbins == 1 & ybins == 1;
            sepoints = xbins == 2 & ybins == 1;
            
            npts = size(obj.points_,2);
            nwpts = sum(nwpoints);
            nepts = sum(nepoints);
            swpts = sum(swpoints);
            septs = sum(sepoints);
            nptsactual = nwpts + nepts + swpts + septs;
            if(npts ~= nptsactual)
                warning("Losing points!")
            end
            
            % Add points to each child
            northwest.points_ = obj.points_(:,nwpoints);
            northeast.points_ = obj.points_(:,nepoints);
            southwest.points_ = obj.points_(:,swpoints);
            southeast.points_ = obj.points_(:,sepoints);
            northwest.Rpoints_ = obj.Rpoints_(:,:,nwpoints);
            northeast.Rpoints_ = obj.Rpoints_(:,:,nepoints);
            southwest.Rpoints_ = obj.Rpoints_(:,:,swpoints);
            southeast.Rpoints_ = obj.Rpoints_(:,:,sepoints);
            
            % Either divide or initialize each child
            %             if(nwpts > obj.nPointsMax_ && (obj.depth_ + 1) < obj.maxDepth_)
            %                 northwest = northwest.Divide();
            %             else
            %                 northwest = northwest.InitFromML();
            %             end
            %             if(nepts > obj.nPointsMax_ && (obj.depth_ + 1) < obj.maxDepth_)
            %                 northeast = northeast.Divide();
            %             else
            %                 northeast = northeast.InitFromML();
            %             end
            %             if(swpts > obj.nPointsMax_ && (obj.depth_ + 1) < obj.maxDepth_)
            %                 southwest = southwest.Divide();
            %             else
            %                 southwest = southwest.InitFromML();
            %             end
            %             if(septs > obj.nPointsMax_ && (obj.depth_ + 1) < obj.maxDepth_)
            %                 southeast = southeast.Divide();
            %             else
            %                 southeast = southeast.InitFromML();
            %             end
            
    
            northwest = northwest.InitFromML();
            northeast = northeast.InitFromML();
            southwest = southwest.InitFromML();
            southeast = southeast.InitFromML();

            
            % Assign children
            obj.northwest_ = northwest;
            obj.northeast_ = northeast;
            obj.southwest_ = southwest;
            obj.southeast_ = southeast;
            
            % Update parameters
            obj.isDivided_ = true;
            obj.isValid_ = true;
            obj.points_ = [];
            obj.Rpoints_ = [];
            
        end
        
        function obj = InitFromML(obj)
            
            % Reset some metrics
            obj.numOutliers_ = 0; % Number of outliers detected
            obj.numReject_ = 0; % Number of measurements rejected
            obj.fitScore_ = 0; % Sum of outlier residuals
            
            % locals
            rmat = obj.points_;
            Rmat = obj.Rpoints_;
            npts = size(rmat,2);
            J = [eye(2); zeros(1,2)];
            
            % Only initialize if npts big enough
            if(npts < 10)
                return
            end
            
            % transform points into cell frame
            rmatcell = rmat;
            rmatcell(1:2,:) = rmatcell(1:2,:) - obj.center_;
            
            % Initialize estimate
            nhat = [0 0 1]';
            d = 0;
            
            % Convergence criteria
            max_iters = 100;
            tol = 1E-6;
            
            % loop
            for ii = 1:max_iters
                
                % find nhatplus
                kappa = nhat(1:2);
                lambda = nhat(3);
                nhatplus = zeros(3);
                nhatplus(1:2,1:2) = eye(2) - kappa*kappa'/(1+lambda);
                nhatplus(3,1:2) = -kappa';
                nhatplus(1:2,3) = kappa;
                nhatplus(3,3) = lambda;
                
                %Find some metrics for each measurement
                b = zeros(3,1);
                H = zeros(3,3);
                for jj = 1:npts
                    cov_eta = nhat'*Rmat(:,:,jj)*nhat;
                    A = rmatcell(:,jj)'*nhat - d;
                    D = [-rmatcell(:,jj)'*nhatplus*J, 1];
                    b = b + D'*A/cov_eta;
                    H = H + D'*D/cov_eta;
                end
                
                % find optimal step
                step = H\b;
                
                % update
                phi = norm(step(1:2));
                if(phi == 0)
                    expphi = [step(1:2); 1];
                else
                    expphi = [sin(phi)*step(1:2)/phi;
                        cos(phi)];
                end
                nhat = nhatplus*expphi;
                d = d + step(3);
                
                % check for convergence
                if(norm(step) < tol)
                    break;
                end
                
            end
            
            if(ii == max_iters)
                warning('No ML convergence detected')
            end
            
            % Assign estimate
            obj.nhat_ = nhat;
            obj.dhat_ = d;
            obj.Phat_ = inv(H);
            obj.isValid_ = true;
            
        end
        
        function PlotTree(obj, fighandle)
            
            % Set figure
            figure(fighandle)
            hold on
            
            % Only continue if tree has valid estimates
            if(~obj.isValid_)
                return;
            end
            
            % Plot this cell or its children
            if(obj.isDivided_)
                obj.northwest_.PlotTree(fighandle);
                obj.northeast_.PlotTree(fighandle);
                obj.southwest_.PlotTree(fighandle);
                obj.southeast_.PlotTree(fighandle);
            else
                
                % Generate x and y points for plotting
                xcenter = obj.center_(1);
                ycenter = obj.center_(2);
                xL = xcenter - obj.halfWidth_;
                xR = xcenter + obj.halfWidth_;
                yU = ycenter + obj.halfHeight_;
                yD = ycenter - obj.halfHeight_;
                xpts = [xL xR xR xL];
                ypts = [yD yD yU yU];
                
                % Get z points
                zpts = zeros(1,4);
                for ii = 1:4
                    zpts(ii) = obj.GetZCoordinate(xpts(ii),ypts(ii));
                end
                
                if(sum(isnan(xpts)) || sum(isnan(ypts)) || sum(isnan(zpts)))
                    warning('NaN Detected!')
                end
                
                % Plot
                patch(xpts,ypts,zpts,rand(1,3),'FaceColor','none','HandleVisibility','off','LineWidth',2);
            end
            
        end
        
        function z = GetZCoordinate(obj,x,y)
            % Inputs in global frame! Outputs in global frame!
            
            % check for division
            if(obj.isDivided_)
                
                % detirmine which quadrant to query
                if(x < obj.center_(1) && y >= obj.center_(2))
                    z = obj.northwest_.GetZCoordinate(x,y);
                elseif (x >= obj.center_(1) && y >= obj.center_(2))
                    z = obj.northeast_.GetZCoordinate(x,y);
                elseif (x < obj.center_(1) && y < obj.center_(2))
                    z = obj.southwest_.GetZCoordinate(x,y);
                elseif (x >= obj.center_(1) && y < obj.center_(2))
                    z = obj.southeast_.GetZCoordinate(x,y);
                else
                    error('Bad coordinate!')
                end
               
            else  
                if(obj.isValid_)
                    %convert inputs to local
                    xlocal = x - obj.center_(1);
                    ylocal = y - obj.center_(2);

                    % locals
                    nhat = obj.nhat_;
                    dhat = obj.dhat_;

                    % z coordinate
                    z = (dhat - nhat(1)*xlocal - nhat(2)*ylocal)/nhat(3);
                else
                    z = NaN;
                end
            end
            
        end
        
        function obj = CountOutliers(obj, rmat, Rmat, bhatmat, m0mat)
            
            % locals
            npts = size(rmat,2);
            nhat = obj.nhat_;
            dhat = obj.dhat_;
            Pbar = obj.Phat_;
            
            % cycle through each point
            for ii = 1:npts
            
                % Transform measurement to cell frame
                rcell = rmat(:,ii);
                rcell(1:2) = rcell(1:2) - obj.center_;
                m0cell = m0mat(:,ii);
                m0cell(1:2) = m0cell(1:2) - obj.center_;
                bhat = bhatmat(:,ii);

                % Predicted measurement
                t_exp = (dhat - nhat'*m0cell)/(nhat'*bhat);
                r_exp = m0cell + t_exp*bhat;

                % Skip if observation angle is too large or expected
                % range negative
                theta = acos(nhat'*bhat)*180/pi;
                if(theta < 110 || t_exp <= 0)
                    continue;
                end

                % Find the measurement jacobian
                H = obj.JacobianEval(bhat(1),bhat(2),bhat(3),dhat, m0cell(1),m0cell(2),m0cell(3),nhat(1),nhat(2),nhat(3));

                % Find innovation covariance
                Pyy = H*Pbar*H' + Rmat(:,:,ii);

                % Find measurement residual
                res = rcell - r_exp;
                
                %eigen decomposition of covariance
                %[V, D] = eig(Pyy);

                %scale eigen vectors by their eigenvalues
                %V(:,1) = V(:,1)*D(1,1);
                %V(:,2) = V(:,2)*D(2,2);

                %solve system
                %b = V\res;

                %sigma value
                %sig = norm(b);
                %obj.fitScore_ = obj.fitScore_ + sig;
                
                % Mahalanobis distance
                mdist = sqrt(res'*(Pyy\res));
                obj.fitScore_ = obj.fitScore_ + mdist;
                
                % update number of outliers
                if(mdist > 6)
%                     disp('Outlier found!')
                    obj.numOutliers_ = obj.numOutliers_ + 1;
                    
%                     disp(obj.center_)
%                     disp(obj.numOutliers_)
                end
            
            end
            
        end

        function mat = QueryMap(obj, gsd, mat, metric)

            % Return if this cell is not valid
            if ~obj.isValid_
                return
            end

            % Check children or return this cell's data
            if obj.isDivided_
                mat = obj.northwest_.QueryMap(gsd,mat,metric);
                mat = obj.northeast_.QueryMap(gsd,mat,metric);
                mat = obj.southwest_.QueryMap(gsd,mat,metric);
                mat = obj.southeast_.QueryMap(gsd,mat,metric);
            else

                % Find the span of this cell in mat
                [matrow, matcol] = size(mat);
                matcenterrow = matrow/2;
                matcentercol = matcol/2;
                cellspanx = obj.halfWidth_/gsd;
                cellspany = obj.halfHeight_/gsd;
                cellcenterx = obj.center_(1)/gsd + matcentercol;
                cellcentery = obj.center_(2)/gsd + matcenterrow;
 
                % Convert to row/col idxs
                leftbound = round(cellcenterx - cellspanx);
                rightbound = round(cellcenterx + cellspanx);
                upbound = round(cellcentery - cellspany); % y-axis reversed
                downbound = round(cellcentery + cellspany); % y-axis reversed

                % Ensure conformity to mat dimensions
                leftbound = max(leftbound, 1);
                rightbound = min(rightbound, matcol);
                upbound = max(upbound, 1);
                downbound = min(downbound, matrow);

                % Query item
                switch metric
                    case 'AvgFitScore'
                        item = obj.fitScore_/size(obj.points_,2);
                    case 'NumPoints'
                        item = size(obj.points_,2);
                    case 'NumReject'
                        item = obj.numReject_;
                    otherwise
                        error("Invalid Query Metric")
                end

                % Write item
                mat(leftbound:rightbound,upbound:downbound) = item;

            end

        end
        
    end
end

