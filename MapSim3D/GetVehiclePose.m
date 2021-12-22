function [x] = GetVehiclePose(t,traj)
%GetVehiclePose returns the true position of the vehicle at time t under
%a specified trajectory. x = [x_pos, y_pos, z_pos, quatMap2Body]'

% number of timesteps
n = length(t);

%initialize output
x = zeros(7,1);

for ii = 1:n
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

            x(1,ii) = a_x + b_x*t(ii) + c_x*t(ii)^2 + d_x*t(ii)^3;
            x(2,ii) = a_y + b_y*t(ii) + c_y*t(ii)^2 + d_y*t(ii)^3;
            x(3,ii) = a_z + b_z*t(ii) + c_z*t(ii)^2 + d_z*t(ii)^3;
            x(4:7,ii) = angle2quat(-pi/2,pi,0,'ZXZ');

        case 'ShortRamp'

            % Initial conditions
            x0 = 10;
            y0 = -10;
            z0 = 100;
            vx0 = -2;
            vy0 = -2;
            vz0 = -2.5;


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

            x(1,ii) = a_x + b_x*t(ii) + c_x*t(ii)^2 + d_x*t(ii)^3;
            x(2,ii) = a_y + b_y*t(ii) + c_y*t(ii)^2 + d_y*t(ii)^3;
            x(3,ii) = a_z + b_z*t(ii) + c_z*t(ii)^2 + d_z*t(ii)^3;
            x(4:7,ii) = angle2quat(-pi/2,pi,0,'ZXZ');

        otherwise
            disp('Error: Invalid Trajectory')
            return
    end

end

end

