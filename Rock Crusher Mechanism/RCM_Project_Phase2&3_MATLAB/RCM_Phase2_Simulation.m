%%% Simulation in m-file %%%

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
th3 = xlsread('Kinematic_Data.xlsx', ['D1:D', num2str(length(t))]);
th4 = xlsread('Kinematic_Data.xlsx', ['G1:G', num2str(length(t))]);
lh7 = xlsread('Kinematic_Data.xlsx', ['J1:J', num2str(length(t))]);
lv7 = xlsread('Kinematic_Data.xlsx', ['M1:M', num2str(length(t))]);
th5 = xlsread('Kinematic_Data.xlsx', ['P1:P', num2str(length(t))]);
th8 = xlsread('Kinematic_Data.xlsx', ['S1:S', num2str(length(t))]);
th6 = xlsread('Kinematic_Data.xlsx', ['V1:V', num2str(length(t))]);

% Position of Joints
O2 = [0; 0];
O4 = [O2_O4_H; -O2_O4_V];
O6 = [-O2_O6_H; -O2_O6_V];
A = [l2*cos(th2)'; l2*sin(th2)'] + O2;
B = [-l4*cos(th4)'; -l4*sin(th4)'] + O4;
C = [l3p*cos(th3+beta)'; l3p*sin(th3+beta)'] + A;
D = [l8*cos(th8)'; l8*sin(th8)'] + O6;
% D = [l5*cos(th5)'; l5*sin(th5)'] + C;
E = [O6_E*cos(th6)'; O6_E*sin(th6)'] + O6;

% Figure of Simulation
flag = 1;
num = 2;    % Number of Periods
count = 0;
j = 1;
gcf = figure('units','normalized','outerposition', [0, 0, 1, 1]);  % Maximize Figure Window
% Options
axis equal;
axis([-100, 40, -100, 30]);
grid on;
box off;
xlabel('$X_{Position}\ (cm)$', 'Interpreter', 'latex');
ylabel('$Y_{Position}\ (cm)$', 'Interpreter', 'latex');
title('$Position\ Diagram$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
hold on;
% Stationary Joints (Joints of Link 1 (Ground))
O2_point = plot(O2(1), O2(2), '.', 'Markersize', 20, 'Color', [0.6350, 0.0780, 0.1840]);
O4_point = plot(O4(1), O4(2), '.', 'Markersize', 20, 'Color', [0.6350, 0.0780, 0.1840]);
O6_point = plot(O6(1), O6(2), '.', 'Markersize', 20, 'Color', [0.6350, 0.0780, 0.1840]);
while true
    for i = 1:length(t)*num
        % Stop running the program when the figure is closed.
        if ishghandle(gcf) == 0
            clc;
            close all;
            flag = 0;
            break;
        end
        
        % Check Array Indices
        if j > length(t)
            count = count + 1;
            if count == num
                break;
            end
            j = j - (num-1)*length(t);
        end
        
        % Joints (Points)
        A_point = plot(A(1, j), A(2, j), '.', 'Markersize', 20, 'Color', [0, 0.4470, 0.7410]);
        B_point = plot(B(1, j), B(2, j), '.', 'Markersize', 20, 'Color', [0.8500, 0.3250, 0.0980]);
        C_point = plot(C(1, j), C(2, j), '.', 'Markersize', 20, 'Color', [0.9290, 0.6940, 0.1250]);
        D_point = plot(D(1, j), D(2, j), '.', 'Markersize', 20, 'Color', [0.4940, 0.1840, 0.5560]);
        E_point = plot(E(1, j), E(2, j), '.', 'Markersize', 20, 'Color', [0.4660, 0.6740, 0.1880]);
        
        % Lines (Links)
        link2 = line([O2(1), A(1, j)], [O2(2), A(2, j)], 'Color', 'black');
        link3_1 = line([A(1, j), B(1, j)], [A(2, j), B(2, j)], 'Color', 'black');
        link3_2 = line([A(1, j), C(1, j)], [A(2, j), C(2, j)], 'Color', 'black');
        link3_3 = line([B(1, j), C(1, j)], [B(2, j), C(2, j)], 'Color', 'black');
        link4 = line([O4(1), B(1, j)], [O4(2), B(2, j)], 'Color', 'black');
        link5 = line([C(1, j), D(1, j)], [C(2, j), D(2, j)], 'Color', 'black');
        link6 = line([O6(1), E(1, j)], [O6(2), E(2, j)], 'Color', 'black');
        link_DE = line([D(1, j), E(1, j)], [D(2, j), E(2, j)], 'Linestyle', '--', 'Color', 'black');
        
        % Time
        period_text = text(-39, 25, ['$Period:\ $', num2str(count+1)],...
                    'Fontsize', 15, 'Interpreter', 'latex');
        time_text = text(-50, 20, ['$Time\ Elapsed:\ $', num2str(t(j)+count*round(t_end, 2)), '$s$'],...
                    'Fontsize', 15, 'Interpreter', 'latex');
        pause(0.08);
        
        % Delete
        if j < length(t)*num 
            delete(A_point);
            delete(B_point);
            delete(C_point);
            delete(D_point);
            delete(E_point);
            delete(link2);
            delete(link3_1);
            delete(link3_2);
            delete(link3_3);
            delete(link4);
            delete(link5);
            delete(link6);
            delete(link_DE);
            delete(period_text);
            delete(time_text);
        end
        
        j = j + 1;
    end
    
    % Stop running the program when time of simulation ends.
    if count == num
        clc;
        close all;
        break
    end
    
    % Stop running the program when the figure is closed.
    if flag == 0
        clc;
        close all;
        break;
    end
end
