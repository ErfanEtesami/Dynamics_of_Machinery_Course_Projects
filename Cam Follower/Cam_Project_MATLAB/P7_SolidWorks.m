%%% Finilized Cam Profile %%%

clc;
clear;
close all;

L = 6;  % Lift (mm)
beta_r = 0.35 * 2 * pi; % Rise Portion (rad)
beta_f = 0.45 * 2 * pi; % Return (Fall) Portion (rad)
dwell = (2*pi - beta_r - beta_f) / 2;   % Half of Dwell Portion (rad)

Offset = 10;    % mm
R0 = 32;    % Radius of Prime Circle (mm)
Rr = 8; % Radius of Roller (mm)

Y0 = sqrt((R0^2) - (Offset^2));

theta = 0:0.1/180*pi:2*pi; % Rotation of Cam (rad)

% Motion Characteristics
y = inf(1, length(theta));
dy = inf(1, length(theta));
ddy = inf(1, length(theta));
dddy = inf(1, length(theta));

% Profile Characteristics
u = inf(1, length(theta));  % X-coordinate of Roller Center
du = inf(1, length(theta)); % X-velocity of Roller Center
ddu = inf(1, length(theta));    % X-Acceleration of Roller Center
v = inf(1, length(theta));  % Y-coordinate of Roller Center
dv = inf(1, length(theta)); % Y-velocity of Roller Center
ddv = inf(1, length(theta));    % Y-Acceleration of Roller Center
dw = inf(1, length(theta));
u_cam = inf(1, length(theta));  % X-coordinate of Cam
v_cam = inf(1, length(theta));  % Y-coordinate of Cam
rho = inf(1, length(theta));    % Curvature of Pitch Curve
rho_cam = inf(1, length(theta));    % Curvature of Cam Profile
Phi = inf(1, length(theta));    % Pressure Angle: cos
phi = inf(1, length(theta));    % Pressure Angle: tan

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
    
    u(i) = (Y0+y(i))*sin(th) + Offset*cos(th);
    v(i) = (Y0+y(i))*cos(th) - Offset*sin(th);
    du(i) = dy(i)*sin(th) + (Y0+y(i)) * cos(th) - Offset*sin(th);
    dv(i) = dy(i)*cos(th) - (Y0+y(i)) * sin(th) - Offset*cos(th);
    ddu(i) = ddy(i)*sin(th) + 2*dy(i)*cos(th) - (Y0+y(i))*sin(th) - Offset*cos(th);
    ddv(i) = ddy(i)*cos(th) - 2*dy(i)*sin(th) - (Y0+y(i))*cos(th) + Offset*sin(th);
    dw(i) = sqrt((du(i)^2) + (dv(i)^2));
    
    u_cam(i) = u(i) + Rr * dv(i)/dw(i);
    v_cam(i) = v(i) - Rr * du(i)/dw(i);
    
    rho(i) = (dw(i)^3)/(du(i)*ddv(i) - dv(i)*ddu(i));
    rho_cam(i) = rho(i) + Rr;
    
    Phi(i) = acosd(-(dv(i)/dw(i))*sin(th) + (du(i)/dw(i))*cos(th));
    phi(i) = atan2d((dy(i)-Offset),(y(i)+Y0));

    i = i + 1;
end

% Write in Excel File
% xlswrite('Cam_Profile.xlsx', u_cam', 'A1:A3601');   % X-coordinate of Cam
% xlswrite('Cam_Profile.xlsx', v_cam', 'B1:B3601');   % Y-coordinate of Cam
% xlswrite('Cam_Profile.xlsx', zeros(36001, 1), 'C1:C3601');  % Z-coordinate of Cam
