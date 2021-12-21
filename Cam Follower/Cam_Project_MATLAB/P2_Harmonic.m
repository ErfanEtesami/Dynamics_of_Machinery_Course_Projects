%%% Harmonic Motion %%%

clc;
clear;
close all;

L = 6;  % Lift (mm)
beta = 0.3 * 2 * pi;    % Rise Portion (rad)
dwell = 0.2 * 2 * pi;   % Half of Dwell Portion (rad)

theta = 0:0.01/180*pi:2*pi; % Rotation of Cam (rad)

y = inf(1, length(theta));  % Displacement
dy = inf(1, length(theta)); % Velocity
ddy = inf(1, length(theta));    % Acceleration
dddy = inf(1, length(theta));   % Jerk

i = 1;
for th = theta
    if th < beta
        y(i) = L/2 * (1 - cos(pi*th/beta));
        dy(i) = pi*L/2/beta * sin(pi*th/beta);
        ddy(i) = (pi^2)*L/2/(beta^2) * cos(pi*th/beta);
        dddy(i) = -(pi^3)*L/2/(beta^3) * sin(pi*th/beta);
    elseif (th >= beta) && (th <= beta + dwell)
        y(i) = L;
        dy(i) = 0;
        ddy(i) = 0;
        dddy(i) = 0;
    elseif (th > beta + dwell) && (th < beta + dwell + beta)
        y(i) = L/2 * (1 + cos(pi*(th-(beta+dwell))/beta));
        dy(i) = -pi*L/2/beta * sin(pi*(th-(beta+dwell))/beta);
        ddy(i) = -(pi^2)*L/2/(beta^2) * cos(pi*(th-(beta+dwell))/beta);
        dddy(i) = (pi^3)*L/2/(beta^3) * sin(pi*(th-(beta+dwell))/beta);
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
title('$Harmonic:\ Displacement\ Diagram$', 'Interpreter', 'latex');
xticks([0, beta/2/pi*180, beta/pi*180, (beta+dwell)/pi*180, (beta+dwell+beta/2)/pi*180,...
        (beta+dwell+beta)/pi*180, 360]);
xticklabels({'0', num2str(beta/2/pi*180), num2str(beta/pi*180), num2str((beta+dwell)/pi*180),...
             num2str((beta+dwell+beta/2)/pi*180), num2str((beta+dwell+beta)/pi*180), '360'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\9.png', '-dpng', '-r500');

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
axis([0, 360, -5, 5]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Velocity\ (mm/s)$', 'Interpreter', 'latex');
title('$Harmonic:\ Velocity\ Diagram$', 'Interpreter', 'latex');
xticks([0, beta/2/pi*180, beta/pi*180, (beta+dwell)/pi*180, (beta+dwell+beta/2)/pi*180,...
        (beta+dwell+beta)/pi*180, 360]);
xticklabels({'0', num2str(beta/2/pi*180), num2str(beta/pi*180), num2str((beta+dwell)/pi*180),...
             num2str((beta+dwell+beta/2)/pi*180), num2str((beta+dwell+beta)/pi*180), '360'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\10.png', '-dpng', '-r500');

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
axis([0, 360, -25/3, 25/3]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Acceleration\ (mm/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Harmonic:\ Acceleration\ Diagram$', 'Interpreter', 'latex');
xticks([0, beta/2/pi*180, beta/pi*180, (beta+dwell)/pi*180, (beta+dwell+beta/2)/pi*180,...
        (beta+dwell+beta)/pi*180, 360]);
xticklabels({'0', num2str(beta/2/pi*180), num2str(beta/pi*180), num2str((beta+dwell)/pi*180),...
             num2str((beta+dwell+beta/2)/pi*180), num2str((beta+dwell+beta)/pi*180), '360'});
yticks([double(-25/3), -6:2:6, double(25/3)]);
yticklabels({'$-25/3$', '-6', '-4', '-2', '0', '2', '4', '6', '$25/3$'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\11.png', '-dpng', '-r500');

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
axis([0, 360, -125/9, 125/9]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Jerk\ (mm/s\textsuperscript{3})$', 'Interpreter', 'latex');
title('$Harmonic:\ Jerk\ Diagram$', 'Interpreter', 'latex');
xticks([0, beta/2/pi*180, beta/pi*180, (beta+dwell)/pi*180, (beta+dwell+beta/2)/pi*180,...
        (beta+dwell+beta)/pi*180, 360]);
xticklabels({'0', num2str(beta/2/pi*180), num2str(beta/pi*180), num2str((beta+dwell)/pi*180),...
             num2str((beta+dwell+beta/2)/pi*180), num2str((beta+dwell+beta)/pi*180), '360'});
yticks([double(-125/9), -10:5:10, double(125/9)]);
yticklabels({'$-125/9$', '-10', '-5', '0', '5', '10', '$125/9$'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\12.png', '-dpng', '-r500');

% Plot: Jerk: Differentiation
b = diff(ddy)./diff(theta);
figure;
plot(theta(1:end-1)/pi*180, b, 'Linewidth', 1.2);
hold on;
xline(beta/2/pi*180, '--');
xline(beta/pi*180);
xline((beta+dwell)/pi*180);
xline((beta+dwell+beta/2)/pi*180, '--');
xline((beta+dwell+beta)/pi*180);
yline(0);
xlim([0, 360]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Jerk\ (mm/s\textsuperscript{3})$', 'Interpreter', 'latex');
title('$Harmonic:\ Jerk\ Diagram: Differentiation$', 'Interpreter', 'latex');
xticks([0, beta/2/pi*180, beta/pi*180, (beta+dwell)/pi*180, (beta+dwell+beta/2)/pi*180,...
        (beta+dwell+beta)/pi*180, 360]);
xticklabels({'0', num2str(beta/2/pi*180), num2str(beta/pi*180), num2str((beta+dwell)/pi*180),...
             num2str((beta+dwell+beta/2)/pi*180), num2str((beta+dwell+beta)/pi*180), '360'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\13.png', '-dpng', '-r500');
