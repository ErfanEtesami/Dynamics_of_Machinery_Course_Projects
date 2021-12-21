%%% Kinematic Plots %%%

clc;
clear;
close all;

% Design Parameters
O2_O4_V = 77;  % Frame: Vertical Distance (cm);
O2_O4_H = 35;  % Frame: Horizontal Distance (cm);
l1 = sqrt((O2_O4_V^2) + (O2_O4_H^2));   % Frame: O2_O4 Distance (cm)
l2 = 18;   % O2_A (cm)
l3 = 96;   % A_B: One Side of Equilateral Triangle (cm)
l3p = 96;  % A_C: One Side of Equilateral Triangle (cm)
b3 = 30; % B_C: Base of Equilateral Triangle (cm)
l4 = 35;   % O4_B (cm)
l5 = 35;   % C_D (cm)
O2_O6_V = 12;  % Frame: Vertical Distance (cm);
O2_O6_H = 60;    % Frame: Vertical Distance (cm);
l7 = sqrt((O2_O6_V^2) + (O2_O6_H^2));   % Frame: O2_O6 Distance (cm)
O6_E = 64; % (cm)
E_D = 8;  % (cm)
l8 = sqrt((O6_E^2) + (E_D^2));  % O6_D (cm)
th1 = 90/180*pi + atan(O2_O4_H / O2_O4_V);  % Frame: O2_O4 Angle (rad)
th7 = 180/180*pi + atan(O2_O6_V / O2_O6_H); % Frame: O2_O6 Angle (rad)
thv7 = 90/180*pi;   % (rad)
thh7 = 0/180*pi;    % (rad)
beta = -18/180*pi; % Apex Angle (cw - rad)
gamma = -atan(E_D / O6_E); % Angle Between l8 and O6_E (cw - rad)

% Motion Characteristic
rpm = 40;   % Maximum Speed of Input Link (O2_A)
range = 180/180*pi; % Allowable Anglular Posiotion of Input Link (O2_A)
a = rpm * (2*pi/60) / (range/2);    % Coefficient of Angular Position
t_end = 2*pi / a;   % Period
t_step = 0.05;  
t = 0:t_step:t_end; % Total Time of One Period

% Kinematic Data
th2 = xlsread('Kinematic_Data.xlsx', ['A1:A', num2str(length(t))]);
omega2 = xlsread('Kinematic_Data.xlsx', ['B1:B', num2str(length(t))]);
alpha2 = xlsread('Kinematic_Data.xlsx', ['C1:C', num2str(length(t))]);
th3 = xlsread('Kinematic_Data.xlsx', ['D1:D', num2str(length(t))]);
omega3 = xlsread('Kinematic_Data.xlsx', ['E1:E', num2str(length(t)-1)]);
alpha3 = xlsread('Kinematic_Data.xlsx', ['F1:F', num2str(length(t)-2)]);
th4 = xlsread('Kinematic_Data.xlsx', ['G1:G', num2str(length(t))]);
omega4 = xlsread('Kinematic_Data.xlsx', ['H1:H', num2str(length(t)-1)]);
alpha4 = xlsread('Kinematic_Data.xlsx', ['I1:I', num2str(length(t)-2)]);
th5 = xlsread('Kinematic_Data.xlsx', ['P1:P', num2str(length(t))]);
omega5 = xlsread('Kinematic_Data.xlsx', ['Q1:Q', num2str(length(t)-1)]);
alpha5 = xlsread('Kinematic_Data.xlsx', ['R1:R', num2str(length(t)-2)]);
th8 = xlsread('Kinematic_Data.xlsx', ['S1:S', num2str(length(t))]);
omega8 = xlsread('Kinematic_Data.xlsx', ['T1:T', num2str(length(t)-1)]);
alpha8 = xlsread('Kinematic_Data.xlsx', ['U1:U', num2str(length(t)-2)]);
th6 = xlsread('Kinematic_Data.xlsx', ['V1:V', num2str(length(t))]);
omega6 = xlsread('Kinematic_Data.xlsx', ['W1:W', num2str(length(t)-1)]);
alpha6 = xlsread('Kinematic_Data.xlsx', ['X1:X', num2str(length(t)-2)]);

% Position of Joints
O2 = [0; 0];
O4 = [O2_O4_H; -O2_O4_V];
O6 = [-O2_O6_H; -O2_O6_V];
A = [l2*cos(th2)'; l2*sin(th2)'] + O2;
B = [-l4*cos(th4)'; -l4*sin(th4)'] + O4;
C = [l3p*cos(th3+beta)'; l3p*sin(th3+beta)'] + A;
C_relative = [l3p*cos(th3+beta)'; l3p*sin(th3+beta)'];
D = [l5*cos(th5)'; l5*sin(th5)'] + C;
D_relative = [l5*cos(th5)'; l5*sin(th5)'];
E = [O6_E*cos(th6)'; O6_E*sin(th6)'] + O6;

% Motion of Joint A
A_x = A(1, :);
A_Vx = diff(A_x) ./ diff(t);
A_Ax = diff(A_Vx) ./ diff(t(1:end-1));
A_y = A(2, :);
A_Vy = diff(A_y) ./ diff(t);
A_Ay = diff(A_Vy) ./ diff(t(1:end-1));
A_V = sqrt(A_Vx.^2 + A_Vy.^2);
A_A = sqrt(A_Ax.^2 + A_Ay.^2);
A_omega = omega2';  % omega2 -> Link2 Angular Velocity
A_omega_check = A_V / l2;
A_alpha = alpha2';    % alpha2 -> Link2 Angular Acceleration

% Motion of Joint B
B_x = B(1, :);
B_Vx = diff(B_x) ./ diff(t);
B_Ax = diff(B_Vx) ./ diff(t(1:end-1));
B_y = B(2, :);
B_Vy = diff(B_y) ./ diff(t);
B_Ay = diff(B_Vy) ./ diff(t(1:end-1));
B_V = sqrt(B_Vx.^2 + B_Vy.^2);
B_A = sqrt(B_Ax.^2 + B_Ay.^2);
B_omega = omega4';  % omega4 -> Link4 Angular Velocity
B_omega_check = B_V / l4;
B_alpha = alpha4';    % alpha4 -> Link4 Angular Acceleration
B_alpha_check = B_A / l4;

% Motion of Joint C
C_x = C(1, :);
C_x_relative = C_relative(1, :);
C_Vx = diff(C_x) ./ diff(t);
C_Vx_relative = diff(C_x_relative) ./ diff(t);
C_Ax = diff(C_Vx) ./ diff(t(1:end-1));
C_Ax_relative = diff(C_Vx_relative) ./ diff(t(1:end-1));
C_y = C(2, :);
C_y_relative = C_relative(2, :);
C_Vy = diff(C_y) ./ diff(t);
C_Vy_relative = diff(C_y_relative) ./ diff(t);
C_Ay = diff(C_Vy) ./ diff(t(1:end-1));
C_Ay_relative = diff(C_Vy_relative) ./ diff(t(1:end-1));
C_V = sqrt(C_Vx.^2 + C_Vy.^2);
C_V_relative = sqrt(C_Vx_relative.^2 + C_Vy_relative.^2);
C_A = sqrt(C_Ax.^2 + C_Ay.^2);
C_A_relative = sqrt(C_Ax_relative.^2 + C_Ay_relative.^2);
C_omega = omega3';  % omega3 -> Link3 Angular Velocity
C_omega_check = C_V_relative / l3p;
C_alpha = alpha3';    % alpha3 -> Link3 Angular Acceleration
C_alpha_check = C_A_relative / l3p;

% Motion of Joint D
D_x = D(1, :);
D_x_relative = D_relative(1, :);
D_Vx = diff(D_x) ./ diff(t);
D_Vx_relative = diff(D_x_relative) ./ diff(t);
D_Ax = diff(D_Vx) ./ diff(t(1:end-1));
D_Ax_relative = diff(D_Vx_relative) ./ diff(t(1:end-1));
D_y = D(2, :);
D_y_relative = D_relative(2, :);
D_Vy = diff(D_y) ./ diff(t);
D_Vy_relative = diff(D_y_relative) ./ diff(t);
D_Ay = diff(D_Vy) ./ diff(t(1:end-1));
D_Ay_relative = diff(D_Vy_relative) ./ diff(t(1:end-1));
D_V = sqrt(D_Vx.^2 + D_Vy.^2);
D_V_relative = sqrt(D_Vx_relative.^2 + D_Vy_relative.^2);
D_A_relative = sqrt(D_Ax_relative.^2 + D_Ay_relative.^2);
D_A = sqrt(D_Ax.^2 + D_Ay.^2);
D_omega = omega5';  % omega5 -> Link5 Angular Velocity
D_omega_check = D_V_relative / l5;
D_alpha = alpha5';  % alpha5 -> Link5 Angular Acceleration
D_alpha_check = D_A_relative / l5;

% Motion of Joint E
E_x = E(1, :);
E_Vx = diff(E_x) ./ diff(t);
E_Ax = diff(E_Vx) ./ diff(t(1:end-1));
E_y = E(2, :);
E_Vy = diff(E_y) ./ diff(t);
E_Ay = diff(E_Vy) ./ diff(t(1:end-1));
E_V = sqrt(E_Vx.^2 + E_Vy.^2);
E_A = sqrt(E_Ax.^2 + E_Ay.^2);
E_omega = omega6';  % omega6 -> Link6 Angular Velocity
E_omega_check = E_V / O6_E;
E_alpha = alpha6';  % alpha6 -> Link6 Angular Acceleration
E_alpha_check = E_A / O6_E;

% Initial Position of Simulation
figure;
axis equal;
axis([-80, 40, -100, 10]);
grid on;
box off;
xlabel('$X_{Position}\ (cm)$', 'Interpreter', 'latex');
ylabel('$Y_{Position}\ (cm)$', 'Interpreter', 'latex');
title('$Position\ Diagram$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
hold on;
O2_point = plot(O2(1), O2(2), '.', 'Markersize', 20, 'Color', [0.6350, 0.0780, 0.1840]);
O4_point = plot(O4(1), O4(2), '.', 'Markersize', 20, 'Color', [0.6350, 0.0780, 0.1840]);
O6_point = plot(O6(1), O6(2), '.', 'Markersize', 20, 'Color', [0.6350, 0.0780, 0.1840]);
A_point = plot(A(1, 1), A(2, 1), '.', 'Markersize', 20, 'Color', [0, 0.4470, 0.7410]);
B_point = plot(B(1, 1), B(2, 1), '.', 'Markersize', 20, 'Color', [0.8500, 0.3250, 0.0980]);
C_point = plot(C(1, 1), C(2, 1), '.', 'Markersize', 20, 'Color', [0.9290, 0.6940, 0.1250]);
D_point = plot(D(1, 1), D(2, 1), '.', 'Markersize', 20, 'Color', [0.4940, 0.1840, 0.5560]);
E_point = plot(E(1, 1), E(2, 1), '.', 'Markersize', 20, 'Color', [0.4660, 0.6740, 0.1880]);
link2 = line([O2(1), A(1, 1)], [O2(2), A(2, 1)], 'Color', 'black');
link3_1 = line([A(1, 1), B(1, 1)], [A(2, 1), B(2, 1)], 'Color', 'black');
link3_2 = line([A(1, 1), C(1, 1)], [A(2, 1), C(2, 1)], 'Color', 'black');
link3_3 = line([B(1, 1), C(1, 1)], [B(2, 1), C(2, 1)], 'Color', 'black');
link4 = line([O4(1), B(1, 1)], [O4(2), B(2, 1)], 'Color', 'black');
link5 = line([C(1, 1), D(1, 1)], [C(2, 1), D(2, 1)], 'Color', 'black');
link6 = line([O6(1), E(1, 1)], [O6(2), E(2, 1)], 'Color', 'black');
link_DE = line([D(1, 1), E(1, 1)], [D(2, 1), E(2, 1)], 'Linestyle', '--', 'Color', 'black');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\1.png', '-dpng', '-r500');

% Plot - Joint A: Position
figure;
plot(A_x, A_y);
grid on;
xlabel('$X\ (cm)$', 'Interpreter', 'latex');
ylabel('$Y\ (cm)$', 'Interpreter', 'latex');
title('$Position\ of\ Joint\ A$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\2.png', '-dpng', '-r500');

% Plot - Joint A: Rotation (Angle)
figure;
plot(t, th2);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Rotation\ (Angle)\ (rad)$', 'Interpreter', 'latex');
title('$Rotation\ of\ Link\ 2$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\3.png', '-dpng', '-r500');

% Plot - Joint A: Speed (Absolute Velocity)
figure;
plot(t(1:end-1), A_V);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Speed\ (cm/s)$', 'Interpreter', 'latex');
title('$Speed\ (Absolute\ Velocity)\ of\ Joint\ A$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\4.png', '-dpng', '-r500');
         
% Plot - Joint A: Angular Velocity and Angular Speed
figure;
plot(t(1:end-1), A_omega);
hold on;
plot(t(1:end-1), A_omega_check, '--');
xlim([0, t_end]);
grid on;
legend('$Angular\ Velocity$', '$Angular\ Speed$', 'Interpreter', 'latex', 'location', 'southwest');
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Angular\ Velocity\ and\ Angular\ Speed\ (rad/s)$', 'Interpreter', 'latex');
title('$Angular\ Velocity\ and\ Angular\ Speed\ of\ Link\ 2$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\5.png', '-dpng', '-r500');

% Plot - Joint A: Absolute Acceleration
figure;
plot(t(1:end-2), A_A);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Acceleration\ (cm/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Absolute\ Acceleration\ of\ Joint\ A$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\6.png', '-dpng', '-r500');

% Plot - Joint A: Angular Acceleration
figure;
plot(t(1:end-2), A_alpha);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Angular\ Acceleration\ (rad/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Angular\ Acceleration\ of\ Link\ 2$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\7.png', '-dpng', '-r500');

% Plot - Joint B: Position
figure;
plot(B_x, B_y);
grid on;
xlabel('$X\ (cm)$', 'Interpreter', 'latex');
ylabel('$Y\ (cm)$', 'Interpreter', 'latex');
title('$Position\ of\ Joint\ B$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\8.png', '-dpng', '-r500');

% Plot - Joint B: Rotation (Angle)
figure;
plot(t, th4);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Rotation\ (Angle)\ (rad)$', 'Interpreter', 'latex');
title('$Rotation\ of\ Link\ 4$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\9.png', '-dpng', '-r500');

% Plot - Joint B: Speed (Absolute Velocity)
figure;
plot(t(1:end-1), B_V);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Speed\ (cm/s)$', 'Interpreter', 'latex');
title('$Speed\ (Absolute\ Velocity)\ of\ Joint\ B$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\10.png', '-dpng', '-r500');
         
% Plot - Joint B: Angular Velocity and Angular Speed
figure;
plot(t(1:end-1), B_omega);
hold on;
plot(t(1:end-1), B_omega_check, '--');
xlim([0, t_end]);
grid on;
legend('$Angular\ Velocity$', '$Angular\ Speed$', 'Interpreter', 'latex', 'location', 'southeast');
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Angular\ Velocity\ and\ Angular\ Speed\ (rad/s)$', 'Interpreter', 'latex');
title('$Angular\ Velocity\ and\ Angular\ Speed\ of\ Link\ 4$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\11.png', '-dpng', '-r500');

% Plot - Joint B: Absolute Acceleration
figure;
plot(t(1:end-2), B_A);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Acceleration\ (cm/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Absolute\ Acceleration\ of\ Joint\ B$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\12.png', '-dpng', '-r500');

% Plot - Joint B: Angular Acceleration
figure;
plot(t(1:end-2), B_alpha);
hold on;
plot(t(1:end-2), B_alpha_check, '--');
xlim([0, t_end]);
grid on;
legend('$Angular\ Acceleration$', '$Angular\ Acceleration\ (Absolute)$',...
       'Interpreter', 'latex', 'location', 'southwest');
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Angular\ Acceleration\ (rad/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Angular\ Acceleration\ of\ Link\ 4$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\13.png', '-dpng', '-r500');

% Plot - Joint C: Position
figure;
plot(C_x, C_y);
grid on;
xlabel('$X\ (cm)$', 'Interpreter', 'latex');
ylabel('$Y\ (cm)$', 'Interpreter', 'latex');
title('$Position\ of\ Joint\ C$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\14.png', '-dpng', '-r500');

% Plot - Joint C: Rotation (Angle)
figure;
plot(t, th3);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Rotation\ (Angle)\ (rad)$', 'Interpreter', 'latex');
title('$Rotation\ of\ Link\ 3$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\15.png', '-dpng', '-r500');

% Plot - Joint C: Speed (Absolute Velocity)
figure;
plot(t(1:end-1), C_V);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Speed\ (cm/s)$', 'Interpreter', 'latex');
title('$Speed\ (Absolute\ Velocity)\ of\ Joint\ C$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\16.png', '-dpng', '-r500');
         
% Plot - Joint C: Angular Velocity and Angular Speed
figure;
plot(t(1:end-1), C_omega);
hold on;
plot(t(1:end-1), C_omega_check, '--');
xlim([0, t_end]);
grid on;
legend('$Angular\ Velocity$', '$Angular\ Speed$', 'Interpreter', 'latex', 'location', 'southwest');
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Angular\ Velocity\ and\ Angular\ Speed\ (rad/s)$', 'Interpreter', 'latex');
title('$Angular\ Velocity\ and\ Angular\ Speed\ of\ Link\ 3$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\17.png', '-dpng', '-r500');

% Plot - Joint C: Absolute Acceleration
figure;
plot(t(1:end-2), C_A);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Acceleration\ (cm/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Absolute\ Acceleration\ of\ Joint\ C$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\18.png', '-dpng', '-r500');

% Plot - Joint C: Angular Acceleration
figure;
plot(t(1:end-2), C_alpha);
hold on;
plot(t(1:end-2), C_alpha_check, '--');
xlim([0, t_end]);
grid on;
legend('$Angular\ Acceleration$', '$Angular\ Acceleration\ (Absolute)$',...
       'Interpreter', 'latex', 'location', 'best');
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Angular\ Acceleration\ (rad/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Angular\ Acceleration\ of\ Link\ 3$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\19.png', '-dpng', '-r500');

% Plot - Joint D: Position
figure;
plot(D_x, D_y);
grid on;
xlabel('$X\ (cm)$', 'Interpreter', 'latex');
ylabel('$Y\ (cm)$', 'Interpreter', 'latex');
title('$Position\ of\ Joint\ D$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\20.png', '-dpng', '-r500');

% Plot - Joint D: Rotation (Angle)
figure;
plot(t, abs(th5));
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Rotation\ (Angle)\ (rad)$', 'Interpreter', 'latex');
title('$Rotation\ of\ Link\ 5$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\21.png', '-dpng', '-r500');

% Plot - Joint D: Speed (Absolute Velocity)
figure;
plot(t(1:end-1), D_V);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Speed\ (cm/s)$', 'Interpreter', 'latex');
title('$Speed\ (Absolute\ Velocity)\ of\ Joint\ D$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\22.png', '-dpng', '-r500');
         
% Plot - Joint D: Angular Velocity and Angular Speed
figure;
plot(t(1:end-1), D_omega);
hold on;
plot(t(1:end-1), D_omega_check, '--');
xlim([0, t_end]);
grid on;
legend('$Angular\ Velocity$', '$Angular\ Speed$', 'Interpreter', 'latex', 'location', 'southeast');
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Angular\ Velocity\ and\ Angular\ Speed\ (rad/s)$', 'Interpreter', 'latex');
title('$Angular\ Velocity\ and\ Angular\ Speed\ of\ Link\ 5$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\23.png', '-dpng', '-r500');

% Plot - Joint D: Absolute Acceleration
figure;
plot(t(1:end-2), D_A);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Acceleration\ (cm/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Absolute\ Acceleration\ of\ Joint\ D$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\24.png', '-dpng', '-r500');

% Plot - Joint D: Angular Acceleration
figure;
plot(t(1:end-2), D_alpha);
hold on;
plot(t(1:end-2), D_alpha_check, '--');
xlim([0, t_end]);
grid on;
legend('$Angular\ Acceleration$', '$Angular\ Accleration\ (Absolute)$',...
       'Interpreter', 'latex', 'location', 'southwest');
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Angular\ Acceleration\ (rad/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Angular\ Acceleration\ of\ Link\ 5$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\25.png', '-dpng', '-r500');

% Plot - Joint E: Position
figure;
plot(E_x, E_y);
grid on;
xlabel('$X\ (cm)$', 'Interpreter', 'latex');
ylabel('$Y\ (cm)$', 'Interpreter', 'latex');
title('$Position\ of\ Joint\ E$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\26.png', '-dpng', '-r500');

% Plot - Joint E: Rotation (Angle)
figure;
plot(t, th6);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Rotation\ (Angle)\ (rad)$', 'Interpreter', 'latex');
title('$Rotation\ of\ Link\ 6$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\27.png', '-dpng', '-r500');

% Plot - Point E: Speed (Absolute Velocity)
figure;
plot(t(1:end-1), E_V);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Speed\ (cm/s)$', 'Interpreter', 'latex');
title('$Speed\ (Absolute\ Velocity)\ of\ Point\ E$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\28.png', '-dpng', '-r500');
         
% Plot - Point E: Angular Velocity and Angular Speed
figure;
plot(t(1:end-1), E_omega);
hold on;
plot(t(1:end-1), E_omega_check, '--');
xlim([0, t_end]);
grid on;
legend('$Angular\ Velocity$', '$Angular\ Speed$', 'Interpreter', 'latex', 'location', 'northeast');
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Angular\ Velocity\ and\ Angular\ Speed\ (rad/s)$', 'Interpreter', 'latex');
title('$Angular\ Velocity\ and\ Angular\ Speed\ of\ Link\ 6$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\29.png', '-dpng', '-r500');

% Plot - Point E: Absolute Acceleration
figure;
plot(t(1:end-2), E_A);
xlim([0, t_end]);
grid on;
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Acceleration\ (cm/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Absolute\ Acceleration\ of\ Point\ E$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\30.png', '-dpng', '-r500');

% Plot - Point E: Angular Acceleration
figure;
plot(t(1:end-2), E_alpha);
hold on;
plot(t(1:end-2), E_alpha_check, '--');
xlim([0, t_end]);
grid on;
legend('$Angular\ Acceleration$', '$Angular\ Accleration\ (Absolute)$',...
       'Interpreter', 'latex', 'location', 'southeast');
xlabel('$Time\ (s)$', 'Interpreter', 'latex');
ylabel('$Angular\ Acceleration\ (rad/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Angular\ Acceleration\ of\ Link\ 6$', 'Interpreter', 'latex');
xticks([0, 0.5, 1, 1.5, 2, t_end]);
xticklabels({'0', '0.5', '1', '1.5', '2', num2str(round(t_end, 2))});
set(gca, 'TickLabelInterpreter', 'latex');
print(gcf, 'G:\University\Mechanics of Machines\Project\Rock Crusher Mechanism\MATLAB\31.png', '-dpng', '-r500');
