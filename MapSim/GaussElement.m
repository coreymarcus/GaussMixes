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
            P = obj.P_mb;
            
            %find midpoint of the line
            avg_s = 0.5*(min_s + max_s);
            avg_x = avg_s/sqrt(1+m_hat^2);
            avg_y = avg_x*m_hat + b_hat;
            
            %approximate the variance in the t direction
            pt_pm = avg_x/sqrt(1+m_hat^2) - (b_hat + m_hat*avg_x - avg_y)*m_hat*(1+m_hat^2)^(-1.5);
            pt_pb = 1/sqrt(1+m_hat^2);
            var_t = [pt_pm, pt_pb]*P*[pt_pm; pt_pb];
            
            %approximate the variance in the s direction
            k_unif = 1.0; %hueristic for scaling uniform variance
            var_s = k_unif*(1/12)*(max_s - min_s)^2;
            
            %now, rotate covariance into the cartesian frame
            P_st = [var_s, 0;
                0, var_t];
            theta = atan(m_hat);
            R_st2xy = [cos(theta), -sin(theta);
                sin(theta), cos(theta)];
            
            %updated variables
            obj.P_xy = R_st2xy*P_st*R_st2xy';
            obj.mu_xy = [avg_x; avg_y];
            
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
            Pxy = Phat*H';
            Pyy = H*Phat*H' + P_y + eye(n)*((Phat(1,1) + xhat(1)).*diag(P_x));
            
            %kalman gain
            K = Pxy/Pyy;
            
            %update
            xhat = xhat + K*(y - H*xhat);
            Phat = (eye(2) - K*H)*Phat*(eye(2) - K*H)' + K*P_y*K';
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
            
            %update number of observations
            obj.Nobs = obj.Nobs + n;
            
        end
        
        %detirmine how likely it was that a measurement came from this
        %gaussian element
        function [p] = GaussEval(obj, x, y, R)
            
            %stack measurement
            z = [x; y];
            
           
            %total covariance
            P = obj.P_xy + R;
            
            %evaluate
            p = 1/sqrt((2*pi)^length(z) * det(P))*exp(-.5*(z - obj.mu_xy)'/P*(x - obj.mu_xy));
            
            %NaN protection
            if(isnan(p))
                p = 0;
            end
            
            
        end
        
        %detirmine how likely it was that a measurement came from this
        %line segment
        function [p] = GaussEvalST(obj, x, y, ~)
            
            %extract local variables
            m_hat = obj.mu_mb(1);
            b_hat = obj.mu_mb(2);
            min_s = obj.s1;
            max_s = obj.s2;
            P = obj.P_mb;
            
            %find midpoint of the line
            avg_s = 0.5*(min_s + max_s);
            avg_x = avg_s/sqrt(1+m_hat^2);
            avg_y = avg_x*m_hat + b_hat;
            
            %approximate the variance in the t direction
            pt_pm = avg_x/sqrt(1+m_hat^2) - (b_hat + m_hat*avg_x - avg_y)*m_hat*(1+m_hat^2)^(-1.5);
            pt_pb = 1/sqrt(1+m_hat^2);
            var_t = [pt_pm, pt_pb]*P*[pt_pm; pt_pb];
            
            %find s and t
            s = x*sqrt(1+m_hat^2);
            t = (b_hat + m_hat*x - y)/sqrt(1+m_hat^2);
            
            %find p_t
            p_t = 1/sqrt(2*pi*var_t)*exp((-.5/var_t)*t^2);
            
            %just approximate p_s
            if((s < max_s+1) && (s > min_s -1))
                p_s = 1;
            else
                p_s = 0;
            end
            
            p = p_s*p_t;
            
            
        end
        
        %detirmine the sigma value of a measurement
        function sig = GaussSig(obj, x, y, R)
            
            %stack measurement
            z = [x; y];
            
            %extract local variables
            mu = obj.mu_xy;
            P = obj.P_xy;
            
            %total covariance
            Q = P + R;
            
            %distance from mean
            res = z - mu;
            
            %eigen decomposition of covariance
            [V, D] = eig(Q);
            
            %scale eigen vectors by their eigenvalues
            V(:,1) = V(:,1)*D(1,1);
            V(:,2) = V(:,2)*D(2,2);
            
            %solve system
            b = V\res;
            
            %sigma value
            sig = norm(b);                       
            
            
        end
        
        function plothandle = PlotElement(obj, handle)
            
            %switch to the current figure
            figure(handle)
            hold on
            
            %sigma value for ellipse
            sigma = 1.5;
            
            % Calculate the eigenvectors and eigenvalues
            [eigenvec, eigenval ] = eig(obj.P_xy);
            
            % Get the index of the largest eigenvector
            [largest_eigenvec_ind_c, ~] = find(eigenval == max(max(eigenval)));
            largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
            
            % Get the largest eigenvalue
            largest_eigenval = max(max(eigenval));
            
            % Get the smallest eigenvector and eigenvalue
            if(largest_eigenvec_ind_c == 1)
                smallest_eigenval = max(eigenval(:,2));
            else
                smallest_eigenval = max(eigenval(:,1));
            end
            
            % Calculate the angle between the x-axis and the largest eigenvector
            angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
            
            % This angle is between -pi and pi.
            % Let's shift it such that the angle is between 0 and 2pi
            if(angle < 0)
                angle = angle + 2*pi;
            end
            
            %create theta
            theta = 0:.01:2*pi;
            
            %create semi-major and minor axes
            a=sigma*sqrt(largest_eigenval);
            b=sigma*sqrt(smallest_eigenval);
            
            % the ellipse in x and y coordinates
            ellipse_x_r  = a*cos( theta );
            ellipse_y_r  = b*sin( theta );
            
            %Define a rotation matrix
            R = [ cos(angle) sin(angle); -sin(angle) cos(angle) ];
            
            %let's rotate the ellipse to some angle phi
            r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
            
            %plot the line segment and elipse
            x1 = obj.s1/sqrt(1+obj.mu_mb(1)^2);
            x2 = obj.s2/sqrt(1+obj.mu_mb(1)^2);
            y1 = obj.mu_mb'*[x1; 1];
            y2 = obj.mu_mb'*[x2; 1];
            plothandle = plot([x1 x2],[y1 y2],'--x',...
                r_ellipse(:,1)+obj.mu_xy(1),r_ellipse(:,2)+obj.mu_xy(2),...
                'LineWidth',2);
            
        end
        
        % Split an element by halving the s domain
        function [obj1, obj2] = SplitElement(obj)
            
            %build new objects
            obj1 = obj;
            obj2 = obj;
            
            %detirmine s-domain
            s1_old = obj.s1;
            s2_old = obj.s2;
            smid = 0.5*(s1_old + s2_old);
            
            obj1.s2 = smid;
            obj2.s1 = smid;
            
            %half the number of observations
            obj1.Nobs = 0.5*obj1.Nobs;
            obj2.Nobs = 0.5*obj2.Nobs;
            
            %update the gaussians of each object
            obj1 = obj1.Line2GaussUpdate();
            obj2 = obj2.Line2GaussUpdate();
            
            
        end
    end
end

