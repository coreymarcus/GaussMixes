function [x] = GetVehiclePose(t,traj)
%GetVehiclePose returns the true position of the vehicle at time t under
%a specified trajectory. x = [x_pos, y_pos, z_pos, quatMap2Body]'

% number of timesteps
n = length(t);

%initialize output
x = zeros(7,n);
switch traj
    case 'Ramp'

        % Initial conditions
        x0 = 100;
        y0 = -100;
        z0 = 2000;
        vx0 = -50;
        vy0 = -50;
        vz0 = -25;


        % Final conditions
        tf = 60;
        xf = 0;
        yf = 0;
        zf = 0;
        vxf = 0;
        vyf = 0;
        vzf = 0;

        % Set up linear system
        a_x = x0;
        b_x = vx0;
        a_y = y0;
        b_y = vy0;
        a_z = z0;
        b_z = vz0;

        A = [tf^2, tf^3;
            2*tf, 3*tf^2];
        B_x = [xf - a_x - b_x*tf; vxf - b_x];
        B_y = [yf - a_y - b_y*tf; vyf - b_y];
        B_z = [zf - a_z - b_z*tf; vzf - b_z];

        z_x = A\B_x;
        z_y = A\B_y;
        z_z = A\B_z;
        c_x = z_x(1);
        d_x = z_x(2);
        c_y = z_y(1);
        d_y = z_y(2);
        c_z = z_z(1);
        d_z = z_z(2);

        for ii = 1:n
            x(1,ii) = a_x + b_x*t(ii) + c_x*t(ii)^2 + d_x*t(ii)^3;
            x(2,ii) = a_y + b_y*t(ii) + c_y*t(ii)^2 + d_y*t(ii)^3;
            x(3,ii) = a_z + b_z*t(ii) + c_z*t(ii)^2 + d_z*t(ii)^3;
            x(4:7,ii) = angle2quat(-pi/2,pi,0,'ZXZ');
        end

    case 'ShortRamp'

        % Initial conditions
        x0 = 10;
        y0 = 10;
        z0 = 150;
        vx0 = -2;
        vy0 = -2;
        vz0 = -2.5;


        % Final conditions
        tf = 60;
        xf = 0;
        yf = 0;
        zf = 15;
        vxf = 0;
        vyf = 0;
        vzf = 0;

        % Set up linear system
        a_x = x0;
        b_x = vx0;
        a_y = y0;
        b_y = vy0;
        a_z = z0;
        b_z = vz0;

        A = [tf^2, tf^3;
            2*tf, 3*tf^2];
        B_x = [xf - a_x - b_x*tf; vxf - b_x];
        B_y = [yf - a_y - b_y*tf; vyf - b_y];
        B_z = [zf - a_z - b_z*tf; vzf - b_z];

        z_x = A\B_x;
        z_y = A\B_y;
        z_z = A\B_z;
        c_x = z_x(1);
        d_x = z_x(2);
        c_y = z_y(1);
        d_y = z_y(2);
        c_z = z_z(1);
        d_z = z_z(2);

        for ii = 1:n
            x(1,ii) = a_x + b_x*t(ii) + c_x*t(ii)^2 + d_x*t(ii)^3;
            x(2,ii) = a_y + b_y*t(ii) + c_y*t(ii)^2 + d_y*t(ii)^3;
            x(3,ii) = a_z + b_z*t(ii) + c_z*t(ii)^2 + d_z*t(ii)^3;
            x(4:7,ii) = angle2quat(-pi/2,pi,0,'ZXZ');
        end

    case 'Scan'

        % Setup scan target points
        Nswitch = 12;
        xmin = -27;
        xmax = 27;
        ymin = -27;
        ymax = 27;

        dx = xmax - xmin;
        dy = ymax - ymin;

        % Get total length
        L = Nswitch*dx + dy;

        % Target points
        Ltargs = zeros(1,2*Nswitch);
        xtargs = zeros(1,2*Nswitch);
        ytargs = zeros(1,2*Nswitch);
        Ltargs(2) = dx;
        xtargs(1) = xmin;
        xtargs(2) = xmax;
        ytargs(1) = ymin;
        ytargs(2) = ymin;
        for ii = 2:Nswitch
            Ltargs(2*ii - 1) = Ltargs(2*ii - 2) + dy/(Nswitch-1);
            Ltargs(2*ii) = Ltargs(2*ii - 1) + dx;

            xtargs(2*ii - 1) = xtargs(2*ii - 2);

            if(mod(ii,2))
                xtargs(2*ii) = xmax;
            else
                xtargs(2*ii) = xmin;
            end

            ytargs(2*ii - 1) = ytargs(2*ii - 2) + dy/(Nswitch-1);
            ytargs(2*ii) = ytargs(2*ii - 1);
        end

        % Target times
        ttargs = t(end).*Ltargs./L;

        % Query points
        xdense = interp1(ttargs,xtargs,t);
        ydense = interp1(ttargs,ytargs,t);
        zdense = zeros(1,n);

        % Determine spacecraft location points
        approachangle = 30;
        range = linspace(490,400,n);
        x(2,:) = -1*range.*cosd(approachangle);
        x(3,:) = range.*sind(approachangle);

        % Find bearing vectors
        b = [xdense; ydense; zdense] - x(1:3,:);

        % Find quaternions
        z = [0 0 1]';
        for ii = 1:n
            % Normalize bearing
            bhat = b(:,ii)/norm(b(:,ii));
            b(:,ii) = bhat;

            % Angle component
            x(4,ii) = 1 + dot(z,bhat);
            x(5:7,ii) = cross(z,bhat);
            x(4:7,ii) = x(4:7,ii)/norm(x(4:7,ii));
        end

    otherwise
        error('Invalid Trajectory')
end

end

