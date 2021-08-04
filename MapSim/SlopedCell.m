classdef SlopedCell
    %SlopedCell Variable resolution DEM cell tree structure
    
    properties
        isDivided_ = false; % Is this cell divided?
        leftChild_; % Child cell if divided
        rightChild_; % Child cell if divided
        isValid_ = false; % Does this cell have a valid estimate?
        x1_; % Left boundary
        x2_; % Right boundary
        mu_mb_; % Slope and intercept estimate
        P_mb_; % Slope and intercept estimate covariance
        points_; % Observations contibuting to this cell
        nPointsMax_; % Number of observations after which division occurs
        depth_; % Depth of this cell ( = 0 for top level, = positive for subsequent levels)
        maxDepth_; % Maximum depth for division
    end
    
    methods
        
        function obj = AddPoints(obj, x, P_x, y, P_y)
            % Adds points to the tree and considers an update
            
            % Check to see if we should add points to this or children?
            if(obj.isDivided_)
                % Bin measurements
                bins = discretize(x,[obj.x1_, 0.5*(obj.x1_ + obj.x2_), obj.x2_]);
                if(sum(bins == 1) > 0)
                    obj.leftChild_ = obj.leftChild_.AddPoints(x(bins == 1), ...
                        P_x(bins == 1, bins == 1),y(bins == 1), ...
                        P_y(bins == 1, bins == 1));
                end
                if(sum(bins == 2) > 0)
                    obj.rightChild_ = obj.rightChild_.AddPoints(x(bins == 2), ...
                        P_x(bins == 2, bins == 2),y(bins == 2), ...
                        P_y(bins == 2, bins == 2));
                end
                
                % No more operations on this cell
                return;
            end
            
            % Ensure measurement vectors are columns
            if(size(x,2) > 1)
                x = x';
                if(size(x,2) > 1)
                    disp('Warning: non-column x')
                end
            end
            if(size(y,2) > 1)
                y = y';
                if(size(y,2) > 1)
                    disp('Warning: non-column y')
                end
            end
            
            % Add points to points_
            obj.points_ = [obj.points_, [x'; y']];
            
            % Decide if we should divide or update
            if((size(obj.points_,2) > obj.nPointsMax_) && obj.depth_ < obj.maxDepth_)
                obj = obj.Divide();
            else
                obj = obj.UpdateKF(x, P_x, y, P_y);
                %                 obj = obj.InitFromTLS();
            end
        end
        
        function obj = UpdateKF(obj, x, P_x, y, P_y)
            
            %extract local variables
            xhat = obj.mu_mb_;
            Phat = obj.P_mb_;
            n = length(y); %number of measurements
            
            % check for valid estimate
            if(~obj.isValid_)
                warning('Attempted to perform KF on invalid cell')
                return
            end
            
            % Ensure measurement vectors are columns
            if(size(x,2) > 1)
                x = x';
                if(size(x,2) > 1)
                    disp('Warning: non-column x')
                end
            end
            if(size(y,2) > 1)
                y = y';
                if(size(y,2) > 1)
                    disp('Warning: non-column y')
                end
            end
            
            %convert points to local frame
            x = x - obj.x1_;
            
            %create H
            H = ones(n,2);
            H(:,1) = x;
            
            %calculate covariances
            Pxy = Phat*H';
            Pyy = H*Phat*H' + P_y + eye(n)*((Phat(1,1) + xhat(1)).*diag(P_x));
            
            %kalman gain
            K = Pxy/Pyy;
            
            %update
            xhat = xhat + K*(y - H*xhat);
            Phat = (eye(2) - K*H)*Phat*(eye(2) - K*H)' + K*P_y*K';
            
            % assign results
            obj.mu_mb_ = xhat;
            obj.P_mb_ = Phat;
            
        end
        
        function obj = Divide(obj)
            
            % Locals
            x1 = obj.x1_;
            x2 = obj.x2_;
            xMid = 0.5*(x1 + x2);
            
            % Create children
            left = SlopedCell();
            right = SlopedCell();
            
            % Children boundaries
            left.x1_ = x1;
            left.x2_ = xMid;
            right.x1_ = xMid;
            right.x2_ = x2;
            
            % Children depth
            left.depth_ = obj.depth_ + 1;
            right.depth_ = obj.depth_ + 1;
            left.maxDepth_ = obj.maxDepth_;
            right.maxDepth_ = obj.maxDepth_;
            
            % Bin measurements
            bins = discretize(obj.points_(1,:),[x1, xMid, x2]);
            
            % Initialize points of cell
            left.points_ = obj.points_(:,bins == 1);
            right.points_ = obj.points_(:,bins == 2);
            left.nPointsMax_ = obj.nPointsMax_;
            right.nPointsMax_ = obj.nPointsMax_;
            
            % Choose to either initialize the cell with TLS or divide again
            if((size(left.points_,2) > left.nPointsMax_) && left.depth_ < left.maxDepth_)
                left = left.Divide();
            else
                left = left.InitFromTLS();
            end
            if((size(right.points_,2) > right.nPointsMax_) && right.depth_ < right.maxDepth_)
                right = right.Divide();
            else
                right = right.InitFromTLS();
            end
            
            
            % Assign children
            obj.isDivided_ = true;
            obj.leftChild_ = left;
            obj.rightChild_ = right;
            
            % Clear out this object
            obj.points_ = [];
            obj.P_mb_ = [];
            obj.mu_mb_ = [];
            
        end
        
        function obj = InitFromTLS(obj)
            
            % Check number of points
            if(size(obj.points_,2) < 3)
                warning('Attempted to initialize with TLS and less than 3 measurements')
                return;
            end
            
            % Transform x points into cell frame
            xEval = obj.points_(1,:) - obj.x1_;
            
            % Use TLS
            % TODO - figure out input variance
            [m_hat, b_hat, P_hat] = TLS(xEval, obj.points_(2,:), 1, 1);
            
            % Initialize estimate
            obj.mu_mb_ = [m_hat, b_hat]';
            obj.P_mb_ = P_hat;
            obj.isValid_ = true;
        end
        
        function PlotTree(obj, fighandle, plotpoints)
            
            % Check if this cell is divided
            if(obj.isDivided_)
                obj.leftChild_.PlotTree(fighandle, plotpoints);
                obj.rightChild_.PlotTree(fighandle, plotpoints);
            else
                figure(fighandle)
                hold on
                
                % Plot data measurements?
                if(plotpoints)
                    scatter(obj.points_(1,:), obj.points_(2,:))
                end
                
                % Plot this if valid
                if(obj.isValid_)
                    xPts = [0 (obj.x2_ - obj.x1_)/2 (obj.x2_ - obj.x1_)];
                    yPts = obj.mu_mb_(1)*xPts + obj.mu_mb_(2);
                    xPts = xPts + obj.x1_;
                    plot(xPts,yPts,'LineWidth',2)
                end
                
                
            end
            
        end
    end
end

