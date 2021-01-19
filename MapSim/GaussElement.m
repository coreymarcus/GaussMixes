classdef GaussElement
    %GAUSSELEMENT class stores all the information for gaussian element
    
    properties
        mu_xy %gaussian mean
        P_xy %gaussian variance
        mu_mb %line mean
        P_mb %line variance
        s1 %lower bound on line segment
        s2 %upper bound on line segment
        Nobs %number of observations
    end
    
    methods
        function obj = GaussElement(Nobs)
            %GAUSSELEMENT Construct an instance of this class
            obj.Nobs = Nobs;
        end
        
        function obj = Line2GaussUpdate(obj)
            %Line2GaussUpdate takes the estimate of the line and uses it to
            %overwrite the estimate of the gaussian
            
            %extract local variables
            m_hat = obj.mu_mb(1);
            b_hat = obj.mu_mb(2);
            min_s = obj.s1;
            max_s = obj.s2;
            P = obj.P_xy;
            
            %find midpoint of the line
            avg_s = 0.5*(min_s + max_s);
            avg_x = avg_s/sqrt(1+m_hat^2);
            avg_y = avg_x*m_hat + b_hat;
            
            %approximate the variance in the t direction
            pt_pm = avg_x/sqrt(1+m_hat^2) - (b_hat + m_hat*avg_x - avg_y)*m_hat*(1+m_hat^2)^(-1.5);
            pt_pb = 1/sqrt(1+m_hat^2);
            var_t = [pt_pm, pt_pb]*P*[pt_pm; pt_pb];
            
            %approximate the variance in the s direction
            k_unif = 0.75; %hueristic for scaling uniform variance
            var_s = k_unif*(1/12)*(max_s - min_s)^2;
            
            %now, rotate covariance into the cartesian frame
            P_st = [var_s, 0;
                0, var_t];
            theta = atan(m_hat_TLS);
            R_st2xy = [cos(theta), -sin(theta);
                sin(theta), cos(theta)];
            
            %updated variables
            obj.P_xy = R_st2xy*P_st*R_st2xy';
            obj.mu_xy = [avg(x); avg(y)];
            
        end
        
        function obj = UpdateLineEstimateKF(obj, x, P_x, y, P_y)
            %Line 2 gauss update performs an update of the line estimate
            %using a modified kalman filter
            
            %extract local variables
            xhat = obj.mu_mb;
            min_s = obj.s1;
            max_s = obj.s2;
            Phat = obj.P_xy;
            n = length(y); %number of measurements
            
            %create H
            H = ones(n,2);
            H(:,1) = x;
            
            %calculate covariances
            Pxy = H*Phat;
            Pyy = H*Phat*H' + P_y + eye(n)*((Phat(1,1) + xhat(1)).*diag(P_x));
            
            %kalman gain
            K = Pxy/Pyy;
            
            %update
            xhat = xhat + K*(y - H*xhat);
            Phat = (eye(2) - K*H)*Phat*(eye(2) - K*H)' + K*R*K';
            obj.mu_mb = xhat;
            obj.P_mb = Phat;
            
            %for all measurements, calculate the s value
            s = zeros(1,n);
            for ii = 1:n
                
                %find s, this should be simplified somehow
                gamma = y(ii) + x(ii)/xhat(1);
                A = [1, -xhat(1);
                    1, 1/xhat(1)];
                inter = A\[xhat(2); gamma];
                xinter = inter(2);
                s(ii) = xinter*sqrt(1+xhat(1)^2);
                
            end
            
            %consider updating s1 and s2
            if(min(s) < min_s)
                obj.s1 = min(s);
            end
            if(max(s) > max_s)
                obj.s2 = max(s);
            end          
            
            
        end
    end
end

