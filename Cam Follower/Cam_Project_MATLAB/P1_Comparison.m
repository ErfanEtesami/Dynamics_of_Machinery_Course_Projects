%%% Comparison between Harmonic, Cycloidal and 8-th Order Polynomial Motion
%%%

clc;
clear;
close all;

L = 6;  % Lift (mm)
beta = 0.30 * 2 * pi;   % Rise Portion (rad)
th = 0:0.01/180*pi:beta; % Rotation of Cam (rad)

a = th / beta;

% Rise: Displacement
y_h_r = @(a) L/2 * (1 - cos(pi*a));
y_c_r = @(a) L * (a - 1/2/pi * sin(2*pi*a));
y_8_r = @(a) L * (6.09755*(a.^3)-20.78040*(a.^5)+26.73155*(a.^6)-...
                  13.60965*(a.^7)+2.56095*(a.^8));
  
% Rise: Velocity             
yd_h_r = @(a) pi*L/2/beta * sin(pi*a);
yd_c_r = @(a) L/beta * (1 - cos(2*pi*a));
yd_8_r = @(a) L/beta * (18.29265*(a.^2)-103.90200*(a.^4)+160.38930*(a.^5)-...
                        95.26755*(a.^6)+20.48760*(a.^7));

% Rise: Acceleration
ydd_h_r = @(a) (pi^2)*L/2/(beta^2) * cos(pi*a);
ydd_c_r = @(a) 2*pi*L/(beta^2) * sin(2*pi*a);
ydd_8_r = @(a) L/(beta^2) * (36.58530*(a)-415.60800*(a.^3)+801.94650*(a.^4)-...
                             571.60530*(a.^5)+143.41320*(a.^6));

% Rise: Jerk
yddd_h_r = @(a) -(pi^3)*L/2/(beta^3) * sin(pi*a);
yddd_c_r = @(a) 4*(pi^2)*L/(beta^3) * cos(2*pi*a);
yddd_8_r = @(a) L/(beta^3) * (36.58530-1246.82400*(a.^2)+3207.78600*(a.^3)-...
                              2858.02650*(a.^4)+860.47920*(a.^5));

% Return (Fall): Displacement
y_h_f = @(a) L/2 * (1 + cos(pi*a));
y_c_f = @(a) L * (1 - a + 1/2/pi * sin(2*pi*a));
y_8_f = @(a) L * (1-2.63415*(a.^2)+2.78055*(a.^5)+3.17060*(a.^6)-...
                  6.87795*(a.^7)+2.56095*(a.^8));
                      
% Return (Fall): Velocity
yd_h_f = @(a) -pi*L/2/beta * sin(pi*a);
yd_c_f = @(a) -L/beta * (1 - cos(2*pi*a));
yd_8_f = @(a) -L/beta * (5.26830*(a)-13.90275*(a.^4)-19.02360*(a.^5)+...
                         48.14565*(a.^6)-20.48760*(a.^7));
   
% Return (Fall): Acceleration
ydd_h_f = @(a) -(pi^2)*L/2/(beta^2) * cos(pi*a);
ydd_c_f = @(a) -2*pi*L/(beta^2) * sin(2*pi*a);
ydd_8_f = @(a) -L/(beta^2) * (5.26830-55.61100*(a.^3)-95.11800*(a.^4)+...
                              288.87390*(a.^5)-143.41320*(a.^6));

% Return (Fall): Jerk
yddd_h_f = @(a) (pi^3)*L/2/(beta^3) * sin(pi*a);
yddd_c_f = @(a) -4*(pi^2)*L/(beta^3) * cos(2*pi*a);
yddd_8_f = @(a) L/(beta^3) * (166.83300*(a.^2)+380.47200*(a.^3)-...
                              1444.36950*(a.^4)+860.47920*(a.^5));                          

% Plot: Rise: Displacement
figure;
fplot(y_h_r, [0, 1]);
hold on;
fplot(y_c_r, [0, 1]);
fplot(y_8_r, [0, 1]);
grid on;
legend('$Harmonic$', '$Cycloidal$', '$8th-order\ Polynomial$',...
       'Interpreter', 'latex', 'Location', 'northwest');
xlabel('$\theta/\beta$', 'Interpreter', 'latex');
ylabel('$Displacement\ (mm)$', 'Interpreter', 'latex');
title('$Rise:\ Displacement\ Compariosn$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\1.png', '-dpng', '-r500');

% Plot: Rise: Velocity
figure;
fplot(yd_h_r, [0, 1]);
hold on;
fplot(yd_c_r, [0, 1]);
fplot(yd_8_r, [0, 1]);
grid on;
legend('$Harmonic$', '$Cycloidal$', '$8th-order\ Polynomial$',...
       'Interpreter', 'latex', 'Units', 'centimeters', 'Position', [7.173, 1.5, 1, 1.3]);
xlabel('$\theta/\beta$', 'Interpreter', 'latex');
ylabel('$Velocity\ (mm/s)$', 'Interpreter', 'latex');
title('$Rise:\ Velocity\ Compariosn$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\2.png', '-dpng', '-r500');

% Plot: Rise: Acceleration
figure;
fplot(ydd_h_r, [0, 1]);
hold on;
fplot(ydd_c_r, [0, 1]);
fplot(ydd_8_r, [0, 1]);
grid on;
yline(0);
legend('$Harmonic$', '$Cycloidal$', '$8th-order\ Polynomial$',...
       'Interpreter', 'latex', 'Location', 'northeast');
xlabel('$\theta/\beta$', 'Interpreter', 'latex');
title('$Rise:\ Acceleration\ Compariosn$', 'Interpreter', 'latex');
ylabel('$Acceleration\ (mm/s\textsuperscript{2})$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\3.png', '-dpng', '-r500');

% Plot: Rise: Jerk
figure;
fplot(yddd_h_r, [0, 1]);
hold on;
fplot(yddd_c_r, [0, 1]);
fplot(yddd_8_r, [0, 1]);
yline(0);
grid on;
legend('$Harmonic$', '$Cycloidal$', '$8th-order\ Polynomial$',...
       'Interpreter', 'latex', 'Units', 'centimeters', 'Position', [7.173, 8.8, 1, 1.3]);
xlabel('$\theta/\beta$', 'Interpreter', 'latex');
ylabel('$Jerk\ (mm/s\textsuperscript{3})$', 'Interpreter', 'latex');
title('$Rise:\ Jerk\ Compariosn$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\4.png', '-dpng', '-r500');

% Plot: Return (Fall): Displacement
figure;
fplot(y_h_f, [0, 1]);
hold on;
fplot(y_c_f, [0, 1]);
fplot(y_8_f, [0, 1]);
grid on;
legend('$Harmonic$', '$Cycloidal$', '$8th-order\ Polynomial$',...
       'Interpreter', 'latex', 'Location', 'northeast');
xlabel('$\theta/\beta$', 'Interpreter', 'latex');
ylabel('$Displacement\ (mm)$', 'Interpreter', 'latex');
title('$Return\ (Fall):\ Displacement\ Compariosn$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\5.png', '-dpng', '-r500');

% Plot: Return (Fall): Velocity
figure;
fplot(yd_h_f, [0, 1]);
hold on;
fplot(yd_c_f, [0, 1]);
fplot(yd_8_f, [0, 1]);
grid on;
legend('$Harmonic$', '$Cycloidal$', '$8th-order\ Polynomial$',...
       'Interpreter', 'latex', 'Units', 'centimeters', 'Position', [7.173, 8.8, 1, 1.3]);
xlabel('$\theta/\beta$', 'Interpreter', 'latex');
ylabel('$Velocity\ (mm/s)$', 'Interpreter', 'latex');
title('$Return\ (Fall):\ Velocity\ Compariosn$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\6.png', '-dpng', '-r500');

% Plot: Return (Fall): Acceleration
figure;
fplot(ydd_h_f, [0, 1]);
hold on;
fplot(ydd_c_f, [0, 1]);
fplot(ydd_8_f, [0, 1]);
grid on;
yline(0);
legend('$Harmonic$', '$Cycloidal$', '$8th-order\ Polynomial$',...
       'Interpreter', 'latex', 'Location', 'northwest');
xlabel('$\theta/\beta$', 'Interpreter', 'latex');
ylabel('$Acceleration\ (mm/s\textsuperscript{2})$', 'Interpreter', 'latex');
title('$Return\ (Fall):\ Acceleration\ Compariosn$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\7.png', '-dpng', '-r500');

% Plot: Return (Fall): Jerk
figure;
fplot(yddd_h_f, [0, 1]);
hold on;
fplot(yddd_c_f, [0, 1]);
fplot(yddd_8_f, [0, 1]);
yline(0);
grid on;
legend('$Harmonic$', '$Cycloidal$', '$8th-order\ Polynomial$',...
       'Interpreter', 'latex', 'Units', 'centimeters', 'Position', [7.173, 1.5, 1, 1.3]);
xlabel('$\theta/\beta$', 'Interpreter', 'latex');
ylabel('$Jerk\ (mm/s\textsuperscript{3})$', 'Interpreter', 'latex');
title('$Return\ (Fall):\ Jerk\ Compariosn$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
% print(gcf, 'G:\University\Mechanics of Machines\Project\Cam & Follower\MATLAB\8.png', '-dpng', '-r500');

% Integrate Accelerations of Rise Portion
i_yh = (yd_h_r(0.5) - yd_h_r(0)) + abs(yd_h_r(1) - yd_h_r(0.5));
i_yc = (yd_c_r(0.5) - yd_c_r(0)) + abs(yd_c_r(1) - yd_c_r(0.5));
r = fzero(ydd_8_r, 0.5);
i_y8 = (yd_8_r(r) - yd_8_r(0)) + abs(yd_8_r(1) - yd_8_r(r));
fprintf('Area Under Acceleration Diagram of Rise: Harmonic = %.4f\n\n', i_yh);
fprintf('Area Under Acceleration Diagram of Rise: Cycloidal = %.4f\n\n', i_yc);
fprintf('Area Under Acceleration Diagram of Rise: 8-th Order Polynomial = %.4f\n\n', i_y8);
