%%% Optimization of Pressure Angle of the Cam %%%

clc;
clear;
close all;

L = 6;  % Lift (mm)
beta_r = 0.35 * 2 * pi; % Rise Portion (rad)
beta_f = 0.45 * 2 * pi; % Return (Fall) Portion (rad)
dwell = (2*pi - beta_r - beta_f) / 2;   % Half of Dwell Portion (rad)

Offset = 10;    % mm
R0 = 15:1:32;    % Radius of Prime Circle (mm)
Rr = 8:0.1:10; % Radius of Roller (mm)

theta = 0:0.01/180*pi:2*pi;

% Motion Characteristics
y = inf(1, length(theta));  % Displacement
dy = inf(1, length(theta)); % Velocity
ddy = inf(1, length(theta));    % Acceleration
dddy = inf(1, length(theta));   % Jerk

i = 1;
for th = theta
    if th < beta_r
        y(i) = L * (th/beta_r - 1/2/pi * sin(2*pi*th/beta_r));
        dy(i) = L/beta_r * (1 - cos(2*pi*th/beta_r));
        ddy(i) = 2*pi*L/(beta_r^2) * sin(2*pi*th/beta_r);
        dddy(i) = 4*(pi^2)*L/(beta_r^3) * cos(2*pi*th/beta_r);
    elseif (th >= beta_r) && (th <= beta_r + dwell)
        y(i) = L;
        dy(i) = 0;
        ddy(i) = 0;
        dddy(i) = 0;
    elseif (th > beta_r + dwell) && (th < beta_r + dwell + beta_f)
        y(i) = L * (1 - (th-(beta_r+dwell))/beta_f + 1/2/pi * sin(2*pi*(th-(beta_r+dwell))/beta_f));
        dy(i) = -L/beta_f * (1 - cos(2*pi*(th-(beta_r+dwell))/beta_f));
        ddy(i) = -2*pi*L/(beta_f^2) * sin(2*pi*(th-(beta_r+dwell))/beta_f);
        dddy(i) = -4*(pi^2)*L/(beta_f^3) * cos(2*pi*(th-(beta_r+dwell))/beta_f);
    elseif th >= beta_r + dwell + beta_f
        y(i) = 0;
        dy(i) = 0;
        ddy(i) = 0;
        dddy(i) = 0;
    end
    
    i = i + 1;
end

% Profile Characteristics
u = inf(length(R0), length(Rr), length(theta)); % X-coordinate of Roller Center
du = inf(length(R0), length(Rr), length(theta));    % X-velocity of Roller Center
ddu = inf(length(R0), length(Rr), length(theta));   % X-Acceleration of Roller Center
v = inf(length(R0), length(Rr), length(theta)); % Y-coordinate of Roller Center
dv = inf(length(R0), length(Rr), length(theta));    % Y-velocity of Roller Center
ddv = inf(length(R0), length(Rr), length(theta));   % Y-Acceleration of Roller Center
dw = inf(length(R0), length(Rr), length(theta));
u_cam = inf(length(R0), length(Rr), length(theta)); % X-coordinate of Cam
v_cam = inf(length(R0), length(Rr), length(theta)); % Y-coordinate of Cam
rho = inf(length(R0), length(Rr), length(theta));   % Curvature of Pitch Curve
rho_cam = inf(length(R0), length(Rr), length(theta));   % Curvature of Cam Profile
Phi = inf(length(R0), length(Rr), length(theta));   % Pressure Angle: cos
phi = inf(length(R0), length(Rr), length(theta));   % Pressure Angle: tan

c1 = 1;
for i = R0
    Y0 = sqrt((i^2) - (Offset^2));
    c2 = 1;
    for j = Rr
        c3 = 1;
        for th = theta
            u(c1, c2, c3) = (Y0+y(c3))*sin(th) + Offset*cos(th);
            v(c1, c2, c3) = (Y0+y(c3))*cos(th) - Offset*sin(th);
            du(c1, c2, c3) = dy(c3)*sin(th) + (Y0+y(c3)) * cos(th) - Offset*sin(th);
            dv(c1, c2, c3) = dy(c3)*cos(th) - (Y0+y(c3)) * sin(th) - Offset*cos(th);
            ddu(c1, c2, c3) = ddy(c3)*sin(th) + 2*dy(c3)*cos(th) - (Y0+y(c3))*sin(th) - Offset*cos(th);
            ddv(c1, c2, c3) = ddy(c3)*cos(th) - 2*dy(c3)*sin(th) - (Y0+y(c3))*cos(th) + Offset*sin(th);
            dw(c1, c2, c3) = sqrt((du(c1, c2, c3)^2) + (dv(c1, c2, c3)^2));
            
            u_cam(c1, c2, c3) = u(c1, c2, c3) + j * dv(c1, c2, c3)/dw(c1, c2, c3);
            v_cam(c1, c2, c3) = v(c1, c2, c3) - j * du(c1, c2, c3)/dw(c1, c2, c3);
    
            rho(c1, c2, c3) = (dw(c1, c2, c3)^3)/(du(c1, c2, c3)*ddv(c1, c2, c3) - dv(c1, c2, c3)*ddu(c1, c2, c3));
            rho_cam(c1, c2, c3) = rho(c1, c2, c3) + j;
    
            Phi(c1, c2, c3) = acosd(-(dv(c1, c2, c3)/dw(c1, c2, c3))*sin(th) + (du(c1, c2, c3)/dw(c1, c2, c3))*cos(th));
            phi(c1, c2, c3) = atan2d((dy(c3)-Offset),(y(c3)+Y0));
            
            c3 = c3 + 1;
        end
        c2 = c2 + 1;
    end
    c1 = c1 + 1;
end

% Find Minimum Amount of Pressure Angle
sum_Phi = inf(c1-1, c2-1, 1);
sum_phi = inf(c1-1, c2-1, 1);
for i = 1:c1-1
    for j = 1:c2-1
        sum_Phi(i, j, 1) = sum(abs(Phi(i, j, :)));
        sum_phi(i, j, 1) = sum(abs(phi(i, j, :)));
    end
end

[min_val, min_loc] = min(sum_Phi(:));
[ii, jj] = ind2sub(size(sum_Phi), min_loc);
[min_val_p, min_loc_p] = min(sum_phi(:));
[ii_p, jj_p] = ind2sub(size(sum_phi), min_loc_p);

fprintf('Radius of Prime Circle: R0 = %0.4f mm\n\n', R0(ii));
fprintf('Radius of Roller: Rr = %0.4f mm\n\n', Rr(jj));
