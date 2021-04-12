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
            k_unif = 1; %hueristic for scaling uniform variance
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
        
        function obj = Gauss2LineUpdate(obj)
            %Gauss2LineUpdate takes the estimate of the gaussian and uses
            %it to overwrite the estimate of the line
            
            %extract local variables
            P = obj.P_xy;
            mu = obj.mu_xy;
            avg_x = mu(1);
            avg_y = mu(2);
            
            %we now need to perform an eigen decomposition on P to
            %detirmine the larger scaled eigenvector
            [V, D] = eig(P);
            v1 = V(:,1)*D(1,1);
            v2 = V(:,2)*D(2,2);
            if(norm(v1) > norm(v2))
                %v1 is the major axis
                m_hat = v1(2)/v1(1);
                theta = atan2(v1(2),v1(1));
                
            elseif(norm(v1) < norm(v2))
                %v2 is the major axis
                m_hat = v2(2)/v2(1);
                theta = atan2(v2(2),v2(1));
                
            else
                %the eigen values are probably the same, so we just take a
                %flat line
                m_hat = 0;
                
            end
            
            %calculate the y intercept
            b_hat = avg_y - m_hat*avg_x;
            
            %rotate the covariance matrix into the s-t domain
            R_xy2st = [cos(theta), -sin(theta);
                sin(theta), cos(theta)]';
            P_st = R_xy2st*P*R_xy2st';
            
            %extract
            var_s = P_st(1,1);
            var_t = P_st(2,2);
            
            %find the s values
            k_unif = 1.0;
            srange = sqrt(12*var_s/k_unif);
            avg_s = avg_x*sqrt(1+m_hat^2);
            min_s = avg_s - 0.5*srange;
            max_s = avg_s + 0.5*srange;
            
            %partial derivatives
            pt_pm = avg_x/sqrt(1+m_hat^2) - (b_hat + m_hat*avg_x - avg_y)*m_hat*(1+m_hat^2)^(-1.5);
            pt_pb = 1/sqrt(1+m_hat^2);
            
            %assume variance of m and b is equivalent
            var_mb = var_t/([pt_pm, pt_pb]*[pt_pm, pt_pb]');
            
            %updated variables
            obj.s1 = min_s;
            obj.s2 = max_s;
            obj.P_mb = var_mb*eye(2);
            obj.mu_mb = [m_hat; b_hat];
            
        end
        
        function obj = UpdateLineEstimateKF(obj, x, P_x, y, P_y)
            %Line 2 gauss update performs an update of the line estimate
            %using a modified kalman filter
            
            %extract local variables
            xhat = obj.mu_mb;
            min_s = obj.s1;
            max_s = obj.s2;
            Phat = obj.P_mb;
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
            s = obj.CalcTransDist(x, y);
            
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
        
        function obj = UpdateGaussDirect(obj, x, P_x, y, P_y)
            %directly updates the gaussian using a sigma point method
            
            %extract local variables
            xhat = obj.mu_xy;
            Phat = obj.P_xy;
            n = length(y); %number of measurements
            sig1 = sqrt(Phat(1,1));
            sig2 = sqrt(Phat(2,2));
            rho = Phat(1,2)/sig1/sig2; %correlation coefficient
            Nmeas = obj.Nobs; %number of measurements already incorporated
            
            sig1_meas = sqrt(P_x);
            sig2_meas = sqrt(P_y);
            rho_meas = 0; %right now uncorrelated error in x and y measurment but this will change!
            
            %cycle through all measurements
            for ii = 1:n
                
                %assemble this measurement as a gaussian
                z = [x(ii), y(ii)]';
                R_meas = zeros(2);
                R_meas(1,1) = P_x;
                R_meas(2,2) = P_y;
                
                %extract sigma points from measurement
                [Xi,w] = GetSigPts(z, R_meas, 0.5);
                
%                 scatter(Xi(1,:),Xi(2,:))
                
                %initialize the updated sigma points
                XiUpdate = zeros(2,5);
                
                % cycle through each sigma point and perform an update
                for jj = 1:5
                    
                    %get the conditional distribution of the measurement
                    %and gaussian conditioned at this sigma point
                    mu_meas = z(2) + rho_meas*sig2_meas*((Xi(1,jj) - z(1))/sig1_meas);
                    var_meas = (1 - rho_meas^2)*sig2_meas^2;
                    
                    mu_gauss = xhat(2) + rho*sig2*((Xi(1,jj) - xhat(1))/sig1);
                    var_gauss = (1 - rho^2)*sig2^2;
                    
                    var_gauss = .1 + var_gauss;
                    
                    %update the sigma point based on the conditionals
%                     XiUpdate(2,jj) = mu_gauss + var_gauss*(var_gauss + var_meas)*(mu_meas - mu_gauss);
                    XiUpdate(2,jj) = mu_gauss + var_gauss*(var_gauss + var_meas)*(Xi(2,jj) - mu_gauss);
                    XiUpdate(1,jj) = Xi(1,jj); %x-component not updated!
                    
                    
                end
                
%                 scatter(XiUpdate(1,:),XiUpdate(2,:),'x')
                
                %find the sample mean and variance of the updated sigma
                %points
                mu_XiUpdate = mean(XiUpdate,2);
                var_XiUpdate = zeros(2);
                diff = XiUpdate - mu_XiUpdate;
                for jj = 1:5
                    var_XiUpdate = var_XiUpdate + w(jj)*diff(:,jj)*diff(:,jj)';
                end
                
                %merge XiUpdate into xhat according to its relative weight
                wXi = 1/(1+Nmeas);
                wGauss = Nmeas/(1+Nmeas);
                
                xhat_update = wXi*mu_XiUpdate + wGauss*xhat;
                Phat = wXi*(var_XiUpdate + mu_XiUpdate*mu_XiUpdate')...
                    + wGauss*(Phat + xhat*xhat') - xhat_update*xhat_update';
                xhat = xhat_update;
                
                Nmeas = Nmeas + 1;
                
                
            end
            
            %update the properties of the object
            obj.mu_xy = xhat;
            obj.P_xy = .00*eye(2) + Phat;
            obj.Nobs = Nmeas;
               
        end
        
        function obj = UpdateGaussNonLinLS(obj, x, P_x, y, P_y)
            %directly updates the gaussian using a non-lin LS optimization
            %scheme
            
            %extract local variables
            xhat = obj.mu_xy;
            Phat = obj.P_xy;
            Nmeas = obj.Nobs;
             
            z = [x, y]';
            R = zeros(2);
            R(1,1) = P_x;
            R(2,2) = P_y;
            
            %form function for optimization
            fun = @(x) NewGaussCost(x,z,xhat,Phat, R, Nmeas);
            
            %form initial guess
            x0 = zeros(5,1);
            x0(1:2) = xhat;
            [V, D] = eig(Phat);
            x0(3) = D(1,1);
            x0(4) = D(2,2);
            x0(5) = V(2,1);
            
            %lower bounds on covariance elements
            LB = [-Inf, -Inf, 0.01, 0.01, -Inf]';
            
            %optimizer
            %  options = optimoptions('lsqnonlin',...
            %   'Algorithm','levenberg-marquardt','Display','off');
            %  xfinal = lsqnonlin(fun,x0,LB,[],options);

            
            options = optimoptions('fmincon','Display','off');
            xfinal = fmincon(fun,x0,[],[],[],[],LB,[],[],options);
            
            fun(x0)
            fun(xfinal)
            
            xhat_new = xfinal(1:2);
            D = zeros(2);
            D(1,1) = xfinal(3);
            D(2,2) = xfinal(4);
            v1 = [1, xfinal(5)]';
            v1 = v1/norm(v1);
            v2 = [-xfinal(5), 1]';
            v2 = v2/norm(v2);
            V = [v1, v2];
            Phat_new = V*D/V;
            
            %update the properties of the object
            obj.mu_xy = xhat_new;
            obj.P_xy = Phat_new;
            obj.Nobs = obj.Nobs + length(x);
               
        end
        
%         function obj = UpdateLineEstimateUKF(obj, x, P_x, y, P_y)
%             %Line 2 gauss update performs an update of the line estimate
%             %using a unscented kalman filter
%             
%             %extract local variables
%             xhat = obj.mu_mb;
%             min_s = obj.s1;
%             max_s = obj.s2;
%             Phat = obj.P_mb;
%             n = length(y); %number of measurements
%             
%             %create H
%             H = ones(n,2);
%             H(:,1) = x;
%             
% %             %calculate covariances
% %             Pxy = Phat*H';
% %             Pyy = H*Phat*H' + P_y + eye(n)*((Phat(1,1) + xhat(1)).*diag(P_x));
% %             
% %             %kalman gain
% %             K = Pxy/Pyy;
% %             
% %             %update
% %             xhat = xhat + K*(y - H*xhat);
% %             Phat = (eye(2) - K*H)*Phat*(eye(2) - K*H)' + K*P_y*K';
% %             obj.mu_mb = xhat;
% %             obj.P_mb = Phat;
% %             
% %             %for all measurements, calculate the s value
% %             s = zeros(1,n);
% %             for ii = 1:n
% %                 
% %                 %find s, this should be simplified somehow
% %                 gamma = y(ii) + x(ii)/xhat(1);
% %                 A = [1, -xhat(1);
% %                     1, 1/xhat(1)];
% %                 inter = A\[xhat(2); gamma];
% %                 xinter = inter(2);
% %                 s(ii) = xinter*sqrt(1+xhat(1)^2);
% %                 
% %             end
% %             
% %             %consider updating s1 and s2
% %             if(min(s) < min_s)
% %                 obj.s1 = min(s);
% %             end
% %             if(max(s) > max_s)
% %                 obj.s2 = max(s);
% %             end
% %             
% %             %update number of observations
% %             obj.Nobs = obj.Nobs + n;
%             
%         end
        
        function obj = UpdateLineEstimateTLS(obj, x, sig2_x, y, sig2_y)
            %Line 2 gauss update performs an update of the line estimate
            %using a modified kalman filter
            
            %extract local variables
            xhat = obj.mu_mb;
            min_s = obj.s1;
            max_s = obj.s2;
            Phat = obj.P_mb;
            n = length(y); %number of measurements
            
            %artificially inflate Phat
            sf = .0;
            Phat = Phat + sf*eye(2);
            
            %perform TLS to get measurements of m and b
            [y_m, y_b, R] = TLS(x, y, sig2_x, sig2_y);
            
            %create H
            H = eye(2);
            
            %calculate covariances
            Pxy = Phat*H';
            Pyy = H*Phat*H' + R;
            
            %kalman gain
            K = Pxy/Pyy;
            
            %update
            xhat = xhat + K*([y_m; y_b] - H*xhat);
            Phat = (eye(2) - K*H)*Phat*(eye(2) - K*H)' + K*R*K';
            obj.mu_mb = xhat;
            obj.P_mb = Phat;
            
            %for all measurements, calculate the s value
            s = obj.CalcTransDist(x, y);
            
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
        
        %determine how likely it was that a measurement came from this
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
        function [p] = GaussEvalST(obj, x, y)
            
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
                p_s = 1/(max_s - min_s);
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
        
        function plothandle = PlotElement(obj, handle, idx)
            
            %switch to the current figure
            figure(handle)
            hold on
            
            %sigma value for ellipse
            sigma = 1;
            
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
            text(x1,y1,string(idx))
            
        end
        
        %calculate the s value for a given point
        function sval = CalcTransDist(obj, x_meas, y_meas)
            
            %extract local variables
            xhat = obj.mu_mb;
            n = length(y_meas); %number of measurements
            
            %for all measurements, calculate the s value
            sval = zeros(1,n);
            for ii = 1:n
                
                %find s, this should be simplified somehow
                gamma = y_meas(ii) + x_meas(ii)/xhat(1);
                A = [1, -xhat(1);
                    1, 1/xhat(1)];
                inter = A\[xhat(2); gamma];
                xinter = inter(2);
                sval(ii) = xinter*sqrt(1+xhat(1)^2);
                
            end
            
            
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
            
            %artificially increase the covariance
            sf = 5;
            obj1.P_mb = sf*obj1.P_mb;
            obj2.P_mb = sf*obj2.P_mb;
            
            
        end
        
    end
end

