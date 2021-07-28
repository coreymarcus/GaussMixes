classdef SlopedCell
    %SlopedCell Variable resolution DEM cell tree structure
    
    properties
        isDivided_ = false; % Is this cell divided?
        leftChild_; % Child cell if divided
        rightChild_; % Child cell if divided
        x1_; % Left boundary
        x2_; % Right boundary
        mu_mb_; % Slope and intercept estimate
        P_mb_; % Slope and intercept estimate covariance
        points_; % Observations contibuting to this cell
        nPointsMax_ = 10; % Number of observations after which division occurs
        maxDepth_ = 5; % Maximum depth for division
    end
    
    methods
        
        function obj = UpdateKF(obj, x, P_x, y, P_y)
            % Updates cell with a kalman filter
            
            
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
            
            % Bin measurements
            bins = discretize(obj.points_(1,:),[x1, xMid, x2]);
            
            % Use TLS to initialize the estimate of each cell
            
            % Assign children
            obj.isDivided_ = true;
            obj.leftChild_ = left;
            obj.rightChild_ = right;
            
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
            [m_hat, b_hat, P_hat] = TLS(xEval, obj.points_(2,:), 1, 1);
            
            % Initialize estimate
            obj.
        end
    end
end

