function [x] = GetVehiclePose(t,traj)
%GetVehiclePose returns the true position of the vehicle at time t under
%a specified trajectory. x = [x_pos, y_pos, theta]'

%initialize output
x = zeros(3,1);

switch traj
    case 'Ramp'
        
        tf = 60;
        x0 = 100; y0 = 1000;
        vx0 = -50; vy0 = -25;
        xf = 0; yf = 0;
        vxf = 0; vyf = 0;
               
        a_x = x0; b_x = vx0;
        a_y = y0; b_y = vy0;
        
        A = [tf^2, tf^3;
            2*tf, 3*tf^2];
        B_x = [xf - a_x - b_x*tf; vxf - b_x];
        B_y = [yf - a_y - b_y*tf; vyf - b_y];
        
        z_x = A\B_x; z_y = A\B_y;
        c_x = z_x(1); d_x = z_x(2);
        c_y = z_y(1); d_y = z_y(2);
        
        x(1) = a_x + b_x*t + c_x*t^2 + d_x*t^3;
        x(2) = a_y + b_y*t + c_y*t^2 + d_y*t^3;
        x(3) = -pi/2;
        
    otherwise
        disp('Error: Invalid Trajectory')
        
end

end

