%%% Cycloidal Motion %%%

clc;
clear;
close all;

L = 6;  % Lift (mm)
beta = 0.3 * 2 * pi;    % Rise Portion (rad)
dwell = 0.2 * 2 * pi;   % Half of Dwell Portion (rad)

theta = 0:0.01/180*pi:2*pi;  % Rotation of Cam (rad)

% Motion Characteristics
y = inf(1, length(theta));  % Displacement
dy = inf(1, length(theta)); % Velocity
ddy = inf(1, length(theta));    % Acceleration
dddy = inf(1, length(theta));   % Jerk

i = 1;
for th = theta
    if th < beta
        y(i) = L * (th/beta - 1/2/pi * sin(2*pi*th/beta));
        dy(i) = L/beta * (1 - cos(2*pi*th/beta));
        ddy(i) = 2*pi*L/(beta^2) * sin(2*pi*th/beta);
        dddy(i) = 4*(pi^2)*L/(beta^3) * cos(2*pi*th/beta);
    elseif (th >= beta) && (th <= beta + dwell)
        y(i) = L;
        dy(i) = 0;
        ddy(i) = 0;
        dddy(i) = 0;
    elseif (th > beta + dwell) && (th < beta + dwell + beta)
        y(i) = L * (1 - (th-(beta+dwell))/beta + 1/2/pi * sin(2*pi*(th-(beta+dwell))/beta));
        dy(i) = -L/beta * (1 - cos(2*pi*(th-(beta+dwell))/beta));
        ddy(i) = -2*pi*L/(beta^2) * sin(2*pi*(th-(beta+dwell))/beta);
        dddy(i) = -4*(pi^2)*L/(beta^3) * cos(2*pi*(th-(beta+dwell))/beta);
    elseif th >= beta + dwell + beta
        y(i) = 0;
        dy(i) = 0;
        ddy(i) = 0;
        dddy(i) = 0;
    end
    
    i = i + 1;
end

% Plot: Displacement
figure;
plot(theta/pi*180, y, 'Linewidth', 1.2);
hold on;
xline(beta/2/pi*180, '--');
xline(beta/pi*180);
xline((beta+dwell)/pi*180);
xline((beta+dwell+beta/2)/pi*180, '--');
xline((beta+dwell+beta)/pi*180);
axis([0, 360, 0, 6]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Displacement\ (mm)$', 'Interpreter', 'latex');
title('$Cycloidal:\ Displacement\ Diagram$', 'Interpreter', 'latex');
xticks([0, beta/2/pi*180, beta/pi*180, (beta+dwell)/pi*180, (beta+dwell+beta/2)/pi*180,...
        (beta+dwell+beta)/pi*180, 360]);
xticklabels({'0', num2str(beta/2/pi*180), num2str(beta/pi*180), num2str((beta+dwell)/pi*180),...
             num2str((beta+dwell+beta/2)/pi*180), num2str((beta+dwell+beta)/pi*180), '360'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\14.png', '-dpng', '-r500');

% Plot: Velocity
figure;
plot(theta/pi*180, dy, 'Linewidth', 1.2);
hold on;
xline(beta/2/pi*180, '--');
xline(beta/pi*180);
xline((beta+dwell)/pi*180);
xline((beta+dwell+beta/2)/pi*180, '--');
xline((beta+dwell+beta)/pi*180);
yline(0);
axis([0, 360, -20/pi, 20/pi]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Velocity\ (mm/s)$', 'Interpreter', 'latex');
title('$Cycloidal:\ Velocity\ Diagram$', 'Interpreter', 'latex');
xticks([0, beta/2/pi*180, beta/pi*180, (beta+dwell)/pi*180, (beta+dwell+beta/2)/pi*180,...
        (beta+dwell+beta)/pi*180, 360]);
xticklabels({'0', num2str(beta/2/pi*180), num2str(beta/pi*180), num2str((beta+dwell)/pi*180),...
             num2str((beta+dwell+beta/2)/pi*180), num2str((beta+dwell+beta)/pi*180), '360'});
yticks([double(-20/pi), -4:2:4, double(20/pi)]);
yticklabels({'$-20/\pi$', '-4', '-2', '0', '2', '4', '$20/\pi$'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\15.png', '-dpng', '-r500');

% Plot: Acceleration
figure;
plot(theta/pi*180, ddy, 'Linewidth', 1.2);
hold on;
xline(beta/2/pi*180, '--');
xline(beta/pi*180);
xline((beta+dwell)/pi*180);
xline((beta+dwell+beta/2)/pi*180, '--');
xline((beta+dwell+beta)/pi*180);
yline(0);
axis([0, 360, -100/3/pi, 100/3/pi]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Acceleration\ (mm/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Cycloidal:\ Acceleration\ Diagram$', 'Interpreter', 'latex');
xticks([0, beta/2/pi*180, beta/pi*180, (beta+dwell)/pi*180, (beta+dwell+beta/2)/pi*180,...
        (beta+dwell+beta)/pi*180, 360]);
xticklabels({'0', num2str(beta/2/pi*180), num2str(beta/pi*180), num2str((beta+dwell)/pi*180),...
             num2str((beta+dwell+beta/2)/pi*180), num2str((beta+dwell+beta)/pi*180), '360'});
yticks([double(-100/3/pi), -8:2:8, double(100/3/pi)]);
yticklabels({'$-100/3\pi$', '-8', '-6', '-4', '-2', '0', '2', '4', '6', '8', '$100/3\pi$'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\16.png', '-dpng', '-r500');

% Plot: Jerk
figure;
plot(theta/pi*180, dddy, 'Linewidth', 1.2);
hold on;
xline(beta/2/pi*180, '--');
xline(beta/pi*180);
xline((beta+dwell)/pi*180);
xline((beta+dwell+beta/2)/pi*180, '--');
xline((beta+dwell+beta)/pi*180);
yline(0);
axis([0, 360, -1000/9/pi, 1000/9/pi]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Jerk\ (mm/s\textsuperscript{3})$', 'Interpreter', 'latex');
title('$Cycloidal:\ Jerk\ Diagram$', 'Interpreter', 'latex');
xticks([0, beta/2/pi*180, beta/pi*180, (beta+dwell)/pi*180, (beta+dwell+beta/2)/pi*180,...
        (beta+dwell+beta)/pi*180, 360]);
xticklabels({'0', num2str(beta/2/pi*180), num2str(beta/pi*180), num2str((beta+dwell)/pi*180),...
             num2str((beta+dwell+beta/2)/pi*180), num2str((beta+dwell+beta)/pi*180), '360'});
yticks([double(-1000/9/pi), -30:10:30, double(1000/9/pi)]);
yticklabels({'$-1000/9\pi$', '-30', '-20', '-10', '0', '10', '20', '30', '$1000/9\pi$'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\17.png', '-dpng', '-r500');
