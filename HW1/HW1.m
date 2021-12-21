clc;
clear;
close all;

a1 = 0.3;   % Frame: OA_OB (m)
a2 = 0.1;   % Crank: OA_A (m)
a3 = 0.2;   % Coupler: A_B (m)
a4 = 0.25;  % Follower: OB_B (m)
a5 = 0.4;   % A_P (m)

t_end = 10;
t_step = 0.05;
t = 0:t_step:t_end;

slope = 5/2;
omega2Max = 5;
t_discontinuity = 2;  % time of discontinuity
t_slope = 0:t_step:t_discontinuity;

omega2 = [slope*t_slope, omega2Max*ones(1,length(t)-length(t_slope))];  % omega_OA_A (Angular Velocity)

th2 = zeros(1, length(t));  % theta_OA_A (Angular Position)
count = 0;
for i = 1:length(t)
    if t(i) <= t_discontinuity
        th2(i) = omega2(i) * t(i) / 2;
        count = count + 1;
    else
        th2(i) = th2(count) + omega2Max*(t(i)-t(count));
    end
end

alpha2 = diff(omega2) ./ diff(t); % alpha_OA_A (Angular Acceleration)

OA = [0; 0];
OB = a1 * [1; 0];
A = [a2*cos(th2); a2*sin(th2)] + OA;
OB_A = sqrt(a1^2 + a2^2 - 2*a1*a2*cos(th2));
beta = asin(a2*sin(th2) ./ OB_A);
gamma = acos((OB_A.^2 + a3^2 - a4^2) ./ (2*a3*OB_A));
B = [a2*cos(th2)+a3*cos(gamma-beta); a2*sin(th2)+a3*sin(gamma-beta)] + OA;
P = [a2*cos(th2)+a5*cos(gamma-beta); a2*sin(th2)+a5*sin(gamma-beta)] + OA;

P_x = P(1, :);
P_Vx = diff(P_x) ./ diff(t);
P_Ax = diff(P_Vx) ./ diff(t(1:end-1));
P_y = P(2, :);
P_Vy = diff(P_y) ./ diff(t);
P_Ay = diff(P_Vy) ./ diff(t(1:end-1));
P_V = sqrt(P_Vx.^2 + P_Vy.^2);
P_A = sqrt(P_Ax.^2 + P_Ay.^2);

delta = acos((OB_A.^2 + a4^2 - a3^2) ./ (2*a4*OB_A));
th4 = -(pi+beta+delta); % theta_OB_B (Angular Position)
omega4 = diff(th4) ./ diff(t);
alpha4 = diff(omega4) ./ diff(t(1:end-1));

scale = 1000;   % m to mm
figure('units','normalized','outerposition',[0, 0, 1, 1]);  % Maximize Figure Window
for i = 1:length(t)
    % Position Diagram
    subplot(3, 2, [1, 2]);
    % Points
    OA_cricle = viscircles(scale*OA', scale*0.003);
    OB_circle = viscircles(scale*OB', scale*0.003);
    A_circle = viscircles(scale*A(:, i)', scale*0.003);
    B_circle = viscircles(scale*B(:, i)', scale*0.003);
    P_circle = viscircles(scale*P(:, i)', scale*0.003);
    % Links
    a1_link = line(scale*[OA(1), OB(1)], scale*[OA(2), OB(2)]);
    a2_link = line(scale*[OA(1), A(1, i)], scale*[OA(2), A(2, i)]);
    a5_link = line(scale*[A(1, i), P(1, i)], scale*[A(2, i), P(2, i)]);
    a4_link = line(scale*[OB(1), B(1, i)], scale*[OB(2), B(2, i)]);
    % Trajectory of Point P n
    hold on;
    plot(scale*P_x(i), scale*P_y(i), 'm.', 'MarkerSize', 3);
    % Options
    grid on;
    axis equal;
    axis(scale*[-0.2, 0.4, -0.2, 0.5]);
    xlabel('$X_{Position}\ (mm)$', 'Interpreter', 'latex');
    ylabel('$Y_{Position}\ (mm)$', 'Interpreter', 'latex');
    title('$Position\ Diagram$', 'Interpreter', 'latex');
    time = text(-scale*0.15, scale*0.45, ['$Time\ Elapsed:\ $', num2str(t(i)), '$s$'], 'Interpreter', 'latex');
    pause(0.005);
    if i < length(t)
        delete(OA_cricle);
        delete(OB_circle);
        delete(A_circle);
        delete(B_circle);
        delete(P_circle);
        delete(a1_link);
        delete(a2_link);
        delete(a5_link);
        delete(a4_link);
        delete(time);
        % Velocity of Point P
        subplot(3, 2, 3);
        plot(t(1:i), P_V(1:i));
        grid on;
        axis([0, 10, 0, 1.5]);
        xlabel('$Time\ (s)$', 'Interpreter', 'latex');
        ylabel('$Velocity\ (m/s)$', 'Interpreter', 'latex');
        title('$Velocity\ of\ Point\ P$', 'Interpreter', 'latex');
        % Angular Velocity of Link 4 (OB_B)
        subplot(3, 2, 4);
        plot(t(1:i), omega4(1:i));
        grid on;
        axis([0, 10, -3, 3]);
        xlabel('$Time\ (s)$', 'Interpreter', 'latex');
        ylabel('$Angular\ Velocity\ (rad/s)$', 'Interpreter', 'latex');
        title('$Angular\ Velocity\ of\ Link\ 4\ (O_B B)$', 'Interpreter', 'latex');
        % Acceleration vectors have to 2 lesser columns than Position vectors
        if i < length(t) - 1    
            % Acceleration of Point P
            subplot(3, 2, 5);
            plot(t(1:i), P_A(1:i));
            grid on;
            axis([0, 10, 0, 8]);
            xlabel('$Time\ (s)$', 'Interpreter', 'latex');
            ylabel('$Acceleration\ (m/s\textsuperscript{2})$', 'Interpreter', 'latex');
            title('$Acceleration\ of\ Point\ P$', 'Interpreter', 'latex');
            % Angular Acceleration of Link 4 (OB_B)
            subplot(3, 2, 6);
            plot(t(1:i), alpha4(1:i));
            grid on;
            axis([0, 10, -20, 20]);
            xlabel('$Time\ (s)$', 'Interpreter', 'latex');
            ylabel('$Angular\ Acceleration\ (rad/s^2)$', 'Interpreter', 'latex');
            title('$Angular\ Acceleration\ of\ Link\ 4\ (O_B B)$', 'Interpreter', 'latex');
        end
    end
end
% print(gcf, 'G:\University\Mechanics of Machines\Zohoor\HWs\MATLAB 1\Final.png', '-dpng', '-r500');

% Check the validity of angular velocity and angular acceleration of OB-B
B_x = B(1, :);
B_Vx = diff(B_x) ./ diff(t);
% B_Ax = diff(B_Vx) ./ diff(t(1:end-1));
B_y = B(2, :);
B_Vy = diff(B_y) ./ diff(t);
% B_Ay = diff(B_Vy) ./ diff(t(1:end-1));
B_V = sqrt(B_Vx.^2 + B_Vy.^2);
% B_A = sqrt(B_Ax.^2 + B_Ay.^2);
omega4_new = B_V / a4;

figure('units','normalized','outerposition',[0, 0, 1, 1]);  % Maximize Figure Window
% True Angular Velocity of Link 4 (OB_B)
subplot(2, 1, 1);
plot(t(1:end-1), omega4);
grid on;
axis([0, 10, -3, 3]);
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Angular\ Velocity\ (rad/s)$', 'Interpreter', 'latex');
title('$True\ Angular\ Velocity\ of\ Link\ 4\ (O_B B)$', 'Interpreter', 'latex');
% Absolute Amount of Angular Velocity of Link 4 (OB_B)
subplot(2, 1, 2);
plot(t(1:end-1), omega4_new);
grid on;
axis([0, 10, 0, 3]);
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Angular\ Velocity\ (rad/s)$', 'Interpreter', 'latex');
title('$abs(Angular\ Velocity\ of\ Link\ 4\ (O_B B))$', 'Interpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Zohoor\HWs\MATLAB 1\Check.png', '-dpng', '-r500');