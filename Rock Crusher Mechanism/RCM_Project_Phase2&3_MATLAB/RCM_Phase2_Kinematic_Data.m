%%% Calculating Kinematic Data %%%

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
th2 = range/2 * (1 -  cos(a*t));    % Angular Position Equation of Input Link (O2_A) Motion
omega2 = diff(th2) ./ diff(t);  % Angular Velocity of Input Link (O2_A)
alpha2 = diff(omega2) ./ diff(t(1:end-1));  % Angular Acceleration of Input Link (O2_A)

% Variables
sol3 = zeros(2, length(th2));   % th3
sol4 = zeros(2, length(th2));   % th4
sol3n = zeros(1, length(th2));  % Accebtable th3
sol4n = zeros(1, length(th2));  % Acceptable th4
sold3 = zeros(1, length(th2));  % omega3
sold4 = zeros(1, length(th2));  % omega4
soldd3 = zeros(1, length(th2)); % alpha3
soldd4 = zeros(1, length(th2)); % alpha4
sollv7 = zeros(1, length(th2)); % lv7
sollh7 = zeros(1, length(th2)); % lh7
soldlv7 = zeros(1, length(th2));    % Velocity of lv7
soldlh7 = zeros(1, length(th2));    % Velocity of lh7
solddlv7 = zeros(1, length(th2));   % Acceleration of lv7
solddlh7 = zeros(1, length(th2));   % Acceleration of lh7  
sol5 = zeros(2, length(th2));   % th5
sol8 = zeros(2, length(th2));   % th8
sol5n = zeros(1, length(th2));  % Acceptable th5
sol8n = zeros(1, length(th2));  % Acceptable th8
sold5 = zeros(1, length(th2));  % omega5
sold8 = zeros(1, length(th2));  % omega8
soldd5 = zeros(1, length(th2)); % alpha5
soldd8 = zeros(1, length(th2)); % alpha8
syms th3 th4;
syms thd3 thd4;
syms thdd3 thdd4;
syms lv7 lh7;
syms lv7d lh7d;
syms lv7dd lh7dd;
syms th5 th8;
syms thd5 thd8;
syms thdd5 thdd8;

% Solving Vector Loop Equations
for i = 1:length(th2)
    vars = [th3, th4];
    eqns = [l2*cos(th2(i)) + l3*cos(th3) + l4*cos(th4) + l1*cos(th1) == 0,...
            l2*sin(th2(i)) + l3*sin(th3) + l4*sin(th4) + l1*sin(th1) == 0];
    [sol3(:, i), sol4(:, i)] = solve(eqns, vars);
    
    sol3n(i) = sol3(1, i);
    sol4n(i) = sol4(1, i);
    
    if i < length(th2)
        vars = [thd3, thd4];
        eqns = [-l2*omega2(i)*sin(th2(i)) - l3*thd3*sin(sol3n(i)) - l4*thd4*sin(sol4n(i)) == 0,...
                l2*omega2(i)*cos(th2(i)) + l3*thd3*cos(sol3n(i)) + l4*thd4*cos(sol4n(i)) == 0];
        [sold3(i), sold4(i)] = solve(eqns, vars);
    end
    
    if i < length(th2) - 1
        vars = [thdd3, thdd4];
        eqns = [-l2*alpha2(i)*sin(th2(i)) - l2*(omega2(i)^2)*cos(th2(i)) - l3*thdd3*sin(sol3n(i)) - l3*(sold3(i)^2)*cos(sol3n(i)) - l4*thdd4*sin(sol4n(i)) - l4*(sold4(i)^2)*cos(sol4n(i)) == 0,...
                l2*alpha2(i)*cos(th2(i)) - l2*(omega2(i)^2)*sin(th2(i)) + l3*thdd3*cos(sol3n(i)) - l3*(sold3(i)^2)*sin(sol3n(i)) + l4*thdd4*cos(sol4n(i)) - l4*(sold4(i)^2)*sin(sol4n(i)) == 0];
        [soldd3(i), soldd4(i)] = solve(eqns, vars);
    end
    
    vars = [lv7, lh7];
    eqns = [-l2*cos(th2(i)) + l7*cos(th7) + lv7*cos(thv7) + lh7*cos(thh7) == 0,...
            -l2*sin(th2(i)) + l7*sin(th7) + lv7*sin(thv7) + lh7*sin(thh7) == 0];
    [sollv7(i), sollh7(i)] = solve(eqns, vars);
    
    if i < length(th2) 
        vars = [lv7d, lh7d];
        eqns = [l2*omega2(i)*sin(th2(i)) + lv7d*cos(thv7) + lh7d*cos(thh7) == 0,...
                -l2*omega2(i)*cos(th2(i)) + lv7d*sin(thv7) + lh7d*sin(thh7) == 0];
        [soldlv7(i), soldlh7(i)] = solve(eqns, vars);
    end
    
    if i < length(th2) - 1
        vars = [lv7dd, lh7dd];
        eqns = [l2*alpha2(i)*sin(th2(i)) - l2*(omega2(i)^2)*cos(th2(i)) + lv7dd*cos(thv7) + lh7dd*cos(thh7) == 0,...
                -l2*alpha2(i)*cos(th2(i)) + l2*(omega2(i)^2)*sin(th2(i)) + lv7dd*sin(thv7) + lh7dd*sin(thh7) == 0];
        [solddlv7(i), solddlh7(i)] = solve(eqns, vars);
    end
    
    vars = [th5, th8];
    eqns = [sollv7(i)*cos(thv7) + sollh7(i)*cos(thh7) + l3p*cos(sol3n(i)+beta) + l5*cos(th5) - l8*cos(th8) == 0,...
            sollv7(i)*sin(thv7) + sollh7(i)*sin(thh7) + l3p*sin(sol3n(i)+beta) + l5*sin(th5) - l8*sin(th8) == 0];
    [sol5(:, i), sol8(:, i)] = solve(eqns, vars);
    
    sol5n(i) = sol5(2, i);
    sol8n(i) = sol8(2, i);
    
    if i < length(th2) 
        vars = [thd5, thd8];
        eqns = [soldlv7(i)*cos(thv7) + soldlh7(i)*cos(thh7) - l3p*sold3(i)*sin(sol3n(i)+beta) - l5*thd5*sin(sol5n(i)) + l8*thd8*sin(sol8n(i)) == 0,...
                soldlv7(i)*sin(thv7) + soldlh7(i)*sin(thh7) + l3p*sold3(i)*cos(sol3n(i)+beta) + l5*thd5*cos(sol5n(i)) - l8*thd8*cos(sol8n(i))];
        [sold5(i), sold8(i)] = solve(eqns, vars);
    end
    
    if i < length(th2) - 1
        vars = [thdd5, thdd8];
        eqns = [solddlv7(i)*cos(thv7) + solddlh7(i)*cos(thh7) - l3p*soldd3(i)*sin(sol3n(i)+beta) - l3p*(sold3(i)^2)*cos(sol3n(i)+beta) - l5*thdd5*sin(sol5n(i)) - l5*(sold5(i)^2)*cos(sol5n(i)) + l8*thdd8*sin(sol8n(i)) + l8*(sold8(i)^2)*cos(sol8n(i)) == 0,...
                solddlv7(i)*sin(thv7) + solddlh7(i)*sin(thh7) + l3p*soldd3(i)*cos(sol3n(i)+beta) - l3p*(sold3(i)^2)*sin(sol3n(i)+beta) + l5*thdd5*cos(sol5n(i)) - l5*(sold5(i)^2)*sin(sol5n(i)) - l8*thdd8*cos(sol8n(i)) + l8*(sold8(i)^2)*sin(sol8n(i)) == 0];
        [soldd5(i), soldd8(i)] = solve(eqns, vars);
    end
end

th6 = sol8n + gamma;
thd6 = diff(th6) ./ diff(t);
thdd6 = diff(thd6) ./ diff(t(1:end-1));

%{
% Write in Excel File
xlswrite('Kinematic_Data.xlsx', th2', ['A1:A', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', omega2', ['B1:B', num2str(length(th2)-1)]);
xlswrite('Kinematic_Data.xlsx', alpha2', ['C1:C', num2str(length(th2)-2)]);
xlswrite('Kinematic_Data.xlsx', sol3', ['D1:D', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sold3', ['E1:E', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', soldd3', ['F1:F', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sol4', ['G1:G', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sold4', ['H1:H', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', soldd4', ['I1:I', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sollh7', ['J1:J', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', soldlh7', ['K1:K', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', solddlh7', ['L1:L', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sollv7', ['M1:M', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', soldlv7', ['N1:N', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', solddlv7', ['O1:O', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sol5n', ['P1:P', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sold5', ['Q1:Q', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', soldd5', ['R1:R', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sol8n', ['S1:S', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sold8', ['T1:T', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', soldd8', ['U1:U', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', th6', ['V1:V', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', thd6', ['W1:W', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', thdd6', ['X1:X', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sol3(2, :)', ['Y1:Y', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sol4(2, :)', ['Z1:Z', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sol5(1, :)', ['AA1:AA', num2str(length(th2))]);
xlswrite('Kinematic_Data.xlsx', sol8(1, :)', ['AB1:AB', num2str(length(th2))]);
%}