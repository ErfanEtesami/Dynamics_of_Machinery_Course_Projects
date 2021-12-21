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

theta = 0:0.01/180*pi:2*pi; % Rotation of Cam (rad)

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

% Plot: Displacement
figure(1);
plot(theta/pi*180, y, 'Linewidth', 1.2);
hold on;
xline(beta_r/2/pi*180, '--');
xline(beta_r/pi*180);
xline((beta_r+dwell)/pi*180);
xline((beta_r+dwell+beta_f/2)/pi*180, '--');
xline((beta_r+dwell+beta_f)/pi*180);
axis([0, 360, 0, 6]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Displacement\ (mm)$', 'Interpreter', 'latex');
title('$Cycloidal\ Displacement\ Diagram\ of\ Follower$', 'Interpreter', 'latex');
xticks([0, beta_r/2/pi*180, beta_r/pi*180, (beta_r+dwell)/pi*180, (beta_r+dwell+beta_f/2)/pi*180,...
        (beta_r+dwell+beta_f)/pi*180, 360]);
xticklabels({'0', num2str(beta_r/2/pi*180), num2str(beta_r/pi*180), num2str((beta_r+dwell)/pi*180),...
             num2str((beta_r+dwell+beta_f/2)/pi*180), num2str((beta_r+dwell+beta_f)/pi*180), '360'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\22.png', '-dpng', '-r500');

% Plot: Velocity
figure(2);
plot(theta/pi*180, dy, 'Linewidth', 1.2);
hold on;
xline(beta_r/2/pi*180, '--');
xline(beta_r/pi*180);
xline((beta_r+dwell)/pi*180);
xline((beta_r+dwell+beta_f/2)/pi*180, '--');
xline((beta_r+dwell+beta_f)/pi*180);
yline(0);
axis([0, 360, -2*L/beta_f, 2*L/beta_r]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Velocity\ (mm/s)$', 'Interpreter', 'latex');
title('$Cycloidal\ Velocity\ Diagram\ of\ Follower$', 'Interpreter', 'latex');
xticks([0, beta_r/2/pi*180, beta_r/pi*180, (beta_r+dwell)/pi*180, (beta_r+dwell+beta_f/2)/pi*180,...
        (beta_r+dwell+beta_f)/pi*180, 360]);
xticklabels({'0', num2str(beta_r/2/pi*180), num2str(beta_r/pi*180), num2str((beta_r+dwell)/pi*180),...
             num2str((beta_r+dwell+beta_f/2)/pi*180), num2str((beta_r+dwell+beta_f)/pi*180), '360'});
yticks([-double(2*L/beta_f), -3:1:4, double(2*L/beta_r)]);
yticklabels({'$-600/(45\pi)$', '-3', '-2', '-1', '0', '1', '2', '3', '4', '$600/(35\pi)$'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\23.png', '-dpng', '-r500');

% Plot: Acceleration
figure(3);
plot(theta/pi*180, ddy, 'Linewidth', 1.2);
hold on;
xline(beta_r/2/pi*180, '--');
xline(beta_r/pi*180);
xline((beta_r+dwell)/pi*180);
xline((beta_r+dwell+beta_f/2)/pi*180, '--');
xline((beta_r+dwell+beta_f)/pi*180);
yline(0);
axis([0, 360, -2*pi*L/(beta_r^2), 2*pi*L/(beta_r^2)]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Acceleration\ (mm/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Cycloidal\ Acceleration\ Diagram\ of\ Follower$', 'Interpreter', 'latex');
xticks([0, beta_r/2/pi*180, beta_r/pi*180, (beta_r+dwell)/pi*180, (beta_r+dwell+beta_f/2)/pi*180,...
        (beta_r+dwell+beta_f)/pi*180, 360]);
xticklabels({'0', num2str(beta_r/2/pi*180), num2str(beta_r/pi*180), num2str((beta_r+dwell)/pi*180),...
             num2str((beta_r+dwell+beta_f/2)/pi*180), num2str((beta_r+dwell+beta_f)/pi*180), '360'});
yticks([double(-2*pi*L/(beta_r^2)), -6:2:6, double(2*pi*L/(beta_r^2))]);
yticklabels({'$-(3\times10\textsuperscript{4})/(35\textsuperscript{2}\pi)$', '-6',...
            '-4', '-2', '0', '2', '4', '6', '$(3\times10\textsuperscript{4})/(35\textsuperscript{2}\pi)$'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\24.png', '-dpng', '-r500');

% Plot: Jerk
figure(4);
plot(theta/pi*180, dddy, 'Linewidth', 1.2);
hold on;
xline(beta_r/2/pi*180, '--');
xline(beta_r/pi*180);
xline((beta_r+dwell)/pi*180);
xline((beta_r+dwell+beta_f/2)/pi*180, '--');
xline((beta_r+dwell+beta_f)/pi*180);
yline(0);
axis([0, 360, -4*(pi^2)*L/(beta_r^3), 4*(pi^2)*L/(beta_r^3)]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Jerk\ (mm/s\textsuperscript{3})$', 'Interpreter', 'latex');
title('$Cycloidal\ Jerk\ Diagram\ of\ Follower$', 'Interpreter', 'latex');
xticks([0, beta_r/2/pi*180, beta_r/pi*180, (beta_r+dwell)/pi*180, (beta_r+dwell+beta_f/2)/pi*180,...
        (beta_r+dwell+beta_f)/pi*180, 360]);
xticklabels({'0', num2str(beta_r/2/pi*180), num2str(beta_r/pi*180), num2str((beta_r+dwell)/pi*180),...
             num2str((beta_r+dwell+beta_f/2)/pi*180), num2str((beta_r+dwell+beta_f)/pi*180), '360'});
yticks([double(-4*(pi^2)*L/(beta_r^3)), -10:10:10, double(4*(pi^2)*L/(beta_r^3))]);
yticklabels({'$-(3\times10\textsuperscript{6})/(35\textsuperscript{3}\pi)$', '-10',...
            '0', '10', '$(3\times10\textsuperscript{6})/(35\textsuperscript{3}\pi)$'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\25.png', '-dpng', '-r500');

% Plot: Cam Profile
figure(5);
plot(u_cam, v_cam);
hold on;
plot(u, v, '--');
xline(0, '--');
yline(0, '--');
legend('$Cam\ Profile$', '$Pitch\ Curve$', 'Interpreter', 'latex');
xlim([-55, 55]);
axis equal;
grid on;
xlabel('$x (mm)$', 'Interpreter', 'latex');
ylabel('$y (mm)$', 'Interpreter', 'latex');
title('$Cam\ Profile$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\26.png', '-dpng', '-r500');

% Plot: Cam Profile and Roller
figure(6);
plot(u_cam, v_cam);
hold on;
plot(u, v, '--');
plot(u(1)+Rr*cos(theta), v(1)+Rr*sin(theta), 'Color', [0.4940, 0.1840, 0.5560]);
plot(u(1), v(1), '.', 'Markersize', 12, 'Color', [0.4660, 0.6740, 0.1880]);
xline(0, '--');
xline(Offset, '--');
yline(0, '--');
legend('$Cam\ Profile$', '$Pitch\ Curve$', '$Roller$', 'Interpreter', 'latex');
ylim([-40, 40]);
axis equal;
grid on;
xlabel('$x (mm)$', 'Interpreter', 'latex');
ylabel('$y (mm)$', 'Interpreter', 'latex');
title('$Cam\ Profile\ and\ Roller$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\27.png', '-dpng', '-r500');

% Animation 1: Check Tangency of Cam Profile and Roller Follower
flag = 1;
gcf = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]); % Maximize Figure Window
ylim([-45, 45]);
axis equal manual;
grid on;
box on;
hold on;
plot(0, 0, 'k.', 'Markersize', 15);
xline(0, '--');
yline(0, '--');
xlabel('$x (mm)$', 'Interpreter', 'latex');
ylabel('$y (mm)$', 'Interpreter', 'latex');
title('$Animation\ of\ Rotating\ Cam\ with\ Reciprocating\ Roller\ Follower\ (Check\ Tangency)$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
while true
    for th = 0:1:360
        % Stop running the program when the figure is closed.
        if ishghandle(gcf) == 0
            clc;
            close all;
            flag = 0;
            break;
        end
        
        x_roller = u*cosd(th) + v*sind(th);
        y_roller = -u*sind(th) + v*cosd(th);
        x_cam = u_cam*cosd(th) + v_cam*sind(th);
        y_cam = -u_cam*sind(th) + v_cam*cosd(th);
        
        cam = plot(x_cam, y_cam, 'Color', [0.8500, 0.3250, 0.0980]);
        cam_contact = plot(x_cam(1), y_cam(1), '.', 'Markersize', 12, 'Color', [0, 0.4470, 0.7410]); 
        roller = plot(x_roller(1)+Rr*cos(theta), y_roller(1)+Rr*sin(theta), 'Color', [0.4940, 0.1840, 0.5560]);
        roller_center = plot(x_roller(1), y_roller(1), '.', 'Markersize', 12, 'Color', [0.4660, 0.6740, 0.1880]);
        common_normal = line([0, x_roller(1)], [0, y_roller(1)], 'LineStyle', '--', 'Color', [0.9290, 0.6940, 0.1250]);
        rotation = text(-4.5, 42, ['$\theta\ =\ $', num2str(th), '$^{\circ}$'],...
                        'Interpreter', 'latex', 'Fontsize', 16);
        
        pause(0.1);
        
        delete(cam);
        delete(cam_contact);
        delete(roller);
        delete(roller_center);
        delete(common_normal);
        delete(rotation);
    end
    
    if th == 360
        clc;
        close all;
        break
    end
    
    if flag == 0
        clc;
        close all;
        break;
    end
end

% Animation 2: Simulation
flag = 1;
gcf = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]); % Maximize Figure Window
ylim([-45, 80]);
axis equal manual;
grid on;
box on;
hold on;
plot(0, 0, 'k.', 'Markersize', 15);
xline(0, '--');
yline(0, '--');
xlabel('$x (mm)$', 'Interpreter', 'latex');
ylabel('$y (mm)$', 'Interpreter', 'latex');
title('$Animation\ of\ Rotating\ Cam\ with\ Reciprocating\ Roller\ Follower\ (Simulation)$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
while true
    for th = 0:1:360
        % Stop running the program when the figure is closed.
        if ishghandle(gcf) == 0
            clc;
            close all;
            flag = 0;
            break;
        end
        
        x_follower = u*cosd(th) + v*sind(th);
        y_follower = -u*sind(th) + v*cosd(th);
        x_cam = u_cam*cosd(th) + v_cam*sind(th);
        y_cam = -u_cam*sind(th) + v_cam*cosd(th);
        
        x_roller = Offset;
        for j = 1:length(y_follower)
            if floor(x_follower(j)) == Offset && y_follower(j) > 0
                break;
            end
        end
        y_roller = y_follower(j);
        
        cam = plot(x_cam, y_cam, 'Color', [0.8500, 0.3250, 0.0980]);
        cam_point = plot(x_cam(1), y_cam(1), '.', 'Markersize', 12, 'Color', [0.3010 0.7450 0.9330]);
        roller = plot(x_roller+Rr*cos(theta), y_roller+Rr*sin(theta), 'Color', [0.4940, 0.1840, 0.5560]);
        roller_center = plot(x_roller, y_roller, '.', 'Markersize', 12, 'Color', [0.4660, 0.6740, 0.1880]);
        roller_point = plot(x_roller+Rr*cosd(th), y_roller+Rr*sind(th), '.', 'MarkerSize', 12, 'Color', [0.6350 0.0780 0.1840]);
        roller_path = line([x_roller, x_roller], [y_roller, 80], 'LineStyle', '--', 'Color', [0.9290, 0.6940, 0.1250]);
        rotation_text = text(-6.3, 77, ['$\theta\ =\ $', num2str(th), '$^{\circ}$'],...
                        'Interpreter', 'latex', 'Fontsize', 16);
        offset_text = text(12, 65, '$Offset$', 'Interpreter', 'latex', 'Fontsize', 12);
        set(offset_text, 'Rotation', 90);
        
        pause(0.1);
        
        delete(cam);
        delete(cam_point);
        delete(roller);
        delete(roller_center);
        delete(roller_point);
        delete(roller_path);
        delete(rotation_text);
        delete(offset_text);
    end
    
    if th == 360
        clc;
        close all;
        break
    end
    
    if flag == 0
        clc;
        close all;
        break;
    end
end
