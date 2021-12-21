%%% Optimization of Acceleration of Cycloidal Displacement Diagram %%%

clc;
clear;
close all;

L = 6:0.1:10;  % Lift (mm)
beta_r = [0.3:0.01:0.35] * 2 * pi;  % Rise Portion (rad)
beta_f = [0.1:0.05:0.45] * 2 * pi;  % Return (Fall) Portion (rad)

theta = 0:0.01/180*pi:2*pi; % Rotation of Cam (rad)

% Motion Characteristics
y = inf(length(L), length(beta_r), length(beta_f), length(theta));  % Displacement
dy = inf(length(L), length(beta_r), length(beta_f), length(theta)); % Velocity
ddy = inf(length(L), length(beta_r), length(beta_f), length(theta));    % Acceleration
dddy = inf(length(L), length(beta_r), length(beta_f), length(theta));   % Jerk

c1 = 1;
for i = L
    c2 = 1;
    for j = beta_r
        c3 = 1;
        for k = beta_f
            c4 = 1;
            dwell = (2*pi - j - k) / 2;
            for th = theta
                if th < j
                    y(c1, c2, c3, c4) = i * (th/j - 1/2/pi * sin(2*pi*th/j));
                    dy(c1, c2, c3, c4) = i/j * (1 - cos(2*pi*th/j));
                    ddy(c1, c2, c3, c4) = 2*pi*i/(j^2) * sin(2*pi*th/j);
                    dddy(c1, c2, c3, c4) = 4*(pi^2)*i/(j^3) * cos(2*pi*th/j);
                elseif (th >= j) && (th <= j + dwell)
                    y(c1, c2, c3, c4) = i;
                    dy(c1, c2, c3, c4) = 0;
                    ddy(c1, c2, c3, c4) = 0;
                    dddy(c1, c2, c3, c4) = 0;
                elseif (th > j + dwell) && (th < j + dwell + k)
                    y(c1, c2, c3, c4) = i * (1 - (th-(j+dwell))/k + 1/2/pi * sin(2*pi*(th-(j+dwell))/k));
                    dy(c1, c2, c3, c4) = -i/k * (1 - cos(2*pi*(th-(j+dwell))/k));
                    ddy(c1, c2, c3, c4) = -2*pi*i/(k^2) * sin(2*pi*(th-(j+dwell))/k);
                    dddy(c1, c2, c3, c4) = -4*(pi^2)*i/(k^3) * cos(2*pi*(th-(j+dwell))/k);
                elseif th >= j + dwell + k
                    y(c1, c2, c3, c4) = 0;
                    dy(c1, c2, c3, c4) = 0;
                    ddy(c1, c2, c3, c4) = 0;
                    dddy(c1, c2, c3, c4) = 0;
                end
                c4 = c4 + 1;
            end
            c3 = c3 + 1;
        end
        c2 = c2 + 1;
    end
    c1 = c1 + 1;
end

% Find Minimum Amount of Accleration
sum_dy = inf(c1-1, c2-1, c3-1, 1);
sum_ddy = inf(c1-1, c2-1, c3-1, 1);
for i = 1:c1-1
    for j = 1:c2-1
        for k = 1:c3-1
            sum_dy(i, j, k, 1) = sum(abs(ddy(i, j, k, :)));
            sum_ddy(i, j, k, 1) = sum(abs(ddy(i, j, k, :)));
        end
    end
end

[min_val, min_loc] = min(sum_ddy(:));
[ii, jj, kk] = ind2sub(size(sum_ddy), min_loc);
[min_val_v, min_loc_v] = min(sum_dy(:));
[ii_v, jj_v, kk_v] = ind2sub(size(sum_dy), min_loc_v);

% Return (Fall): 0.45, 0.35 = Rise, 0.30 Portion
y_45 = reshape(y(ii, jj, kk, :), [length(theta), 1, 1, 1]);
y_35 = reshape(y(ii, jj, 6, :), [length(theta), 1, 1, 1]);
y_30 = reshape(y(ii, jj, 5, :), [length(theta), 1, 1, 1]);
dy_45 = reshape(dy(ii, jj, kk, :), [length(theta), 1, 1, 1]);
dy_35 = reshape(dy(ii, jj, 6, :), [length(theta), 1, 1, 1]);
dy_30 = reshape(dy(ii, jj, 5, :), [length(theta), 1, 1, 1]);
ddy_45 = reshape(ddy(ii, jj, kk, :), [length(theta), 1, 1, 1]);
ddy_35 = reshape(ddy(ii, jj, 6, :), [length(theta), 1, 1, 1]);
ddy_30 = reshape(ddy(ii, jj, 5, :), [length(theta), 1, 1, 1]);
dddy_45 = reshape(dddy(ii, jj, kk, :), [length(theta), 1, 1, 1]);
dddy_35 = reshape(dddy(ii, jj, 6, :), [length(theta), 1, 1, 1]);
dddy_30 = reshape(dddy(ii, jj, 5, :), [length(theta), 1, 1, 1]);

% Plot: Displacement
figure;
plot(theta/pi*180, y_45);
hold on;
plot(theta/pi*180, y_35);
plot(theta/pi*180, y_30);
legend('$\beta_f = 0.45$', '$\beta_f = 0.35$', '$\beta_f = 0.30$',...
       'Interpreter', 'latex', 'Location', 'northeast');
axis([0, 360, 0, 6]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Displacement\ (mm)$', 'Interpreter', 'latex');
title('$Cycloidal:\ Displacement\ Diagrams$', 'Interpreter', 'latex');
xticks([0:40:360]);
xticklabels({'0', '40', '80', '120', '160', '200', '240', '280', '320', '360'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\18.png', '-dpng', '-r500');

% Plot: Velocity
figure;
plot(theta/pi*180, dy_45);
hold on;
plot(theta/pi*180, dy_35);
plot(theta/pi*180, dy_30);
legend('$\beta_f = 0.45$', '$\beta_f = 0.35$', '$\beta_f = 0.30$',...
       'Interpreter', 'latex', 'Location', 'northeast');
axis([0, 360, -7, 7]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Velocity\ (mm/s)$', 'Interpreter', 'latex');
title('$Cycloidal:\ Velocity\ Diagrams$', 'Interpreter', 'latex');
xticks([0:40:360]);
xticklabels({'0', '40', '80', '120', '160', '200', '240', '280', '320', '360'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\19.png', '-dpng', '-r500');

% Plot: Acceleration
figure;
plot(theta/pi*180, ddy_45);
hold on;
plot(theta/pi*180, ddy_35);
plot(theta/pi*180, ddy_30);
legend('$\beta_f = 0.45$', '$\beta_f = 0.35$', '$\beta_f = 0.30$',...
       'Interpreter', 'latex', 'Units', 'centimeters', 'Position', [7.173, 8.8, 1, 1.3]);
axis([0, 360, -11, 11]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Acceleration\ (mm/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Cycloidal:\ Acceleration\ Diagrams$', 'Interpreter', 'latex');
xticks([0:40:360]);
xticklabels({'0', '40', '80', '120', '160', '200', '240', '280', '320', '360'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\20.png', '-dpng', '-r500');

% Plot: Jerk
figure;
plot(theta/pi*180, dddy_45);
hold on;
plot(theta/pi*180, dddy_35);
plot(theta/pi*180, dddy_30);
legend('$\beta_f = 0.45$', '$\beta_f = 0.35$', '$\beta_f = 0.30$',...
       'Interpreter', 'latex', 'Location', 'northeast');
axis([0, 360, -40, 40]);
grid on;
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$Jerk\ (mm/s\textsuperscript{3})$', 'Interpreter', 'latex');
title('$Cycloidal:\ Jerk\ Diagrams$', 'Interpreter', 'latex');
xticks([0:40:360]);
xticklabels({'0', '40', '80', '120', '160', '200', '240', '280', '320', '360'});
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\21.png', '-dpng', '-r500');

fprintf('Lift: L = %0.4f mm\n\n', L(ii));
fprintf('Rise Portion: Beta_r = %0.4f of cycle\n\n', beta_r(jj)/2/pi);
fprintf('Return (Fall) Portion: Beta_f = %0.4f of cycle\n\n', beta_f(kk)/2/pi);
