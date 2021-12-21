clc;
clear;
close all;

l1 = 8.458;
l2 = 1.8;
l3 = 9.6;
l3p = 9.6;
l4 = 3.5;
l5 = 3.5;
l7 = 6.119;
l8 = 6.450;

th1 = 114.445 / 180 * pi;
th7 = 191.310 / 180 * pi;
thv7 = 90 / 180 * pi;
thh7 = 0 / 180 * pi;
beta = -18 / 180 * pi;
gamma = -3.576 / 180 * pi;

th2 = 140 / 180 * pi;
thd2 = -1;
thdd2 = -1;

%% Vector Loop #1
syms th3 th4;

vars = [th3, th4];

eqns = [l2*cos(th2) + l3*cos(th3) + l4*cos(th4) + l1*cos(th1) == 0,...
        l2*sin(th2) + l3*sin(th3) + l4*sin(th4) + l1*sin(th1) == 0];
[sol3, sol4] = solve(eqns, vars);

% disp(double(sol3 * 180 / pi));
% disp(double(sol4 * 180 / pi));

for i = 1:2
    if double(sol4(i) * 180 / pi) >= 10
        sol3n = double(sol3(i));
        sol4n = double(sol4(i));
    end
end

% disp(double(sol3n * 180 / pi));
% disp(double(sol4n * 180 / pi));
% disp(double(sol3n));
% disp(double(sol4n));
% disp('--------------------------');

syms thd3 thd4;

vars = [thd3, thd4];

eqns = [-l2*thd2*sin(th2) - l3*thd3*sin(sol3n) - l4*thd4*sin(sol4n) == 0,...
        l2*thd2*cos(th2) + l3*thd3*cos(sol3n) + l4*thd4*cos(sol4n) == 0];
[sold3, sold4] = solve(eqns, vars);

% disp(double(sold3));
% disp(double(sold4));
% disp('--------------------------');

syms thdd3 thdd4;

vars = [thdd3, thdd4];

eqns = [-l2*thdd2*sin(th2) - l2*(thd2^2)*cos(th2) - l3*thdd3*sin(sol3n) - l3*(sold3^2)*cos(sol3n) - l4*thdd4*sin(sol4n) - l4*(sold4^2)*cos(sol4n) == 0,...
        l2*thdd2*cos(th2) - l2*(thd2^2)*sin(th2) + l3*thdd3*cos(sol3n) - l3*(sold3^2)*sin(sol3n) + l4*thdd4*cos(sol4n) - l4*(sold4^2)*sin(sol4n) == 0];
[soldd3, soldd4] = solve(eqns, vars);

% disp(double(soldd3));
% disp(double(soldd4));
% disp('--------------------------');

%% Vector Loop #2
syms lv7 lh7;

vars = [lv7, lh7];

eqns = [-l2*cos(th2) + l7*cos(th7) + lv7*cos(thv7) + lh7*cos(thh7) == 0,...
        -l2*sin(th2) + l7*sin(th7) + lv7*sin(thv7) + lh7*sin(thh7) == 0];
[sollv7, sollh7] = solve(eqns, vars);

% disp(double(sollv7));
% disp(double(sollh7));
% disp('--------------------------');

syms lv7d lh7d;

vars = [lv7d, lh7d];

eqns = [l2*thd2*sin(th2) + lv7d*cos(thv7) + lh7d*cos(thh7) == 0,...
        -l2*thd2*cos(th2) + lv7d*sin(thv7) + lh7d*sin(thh7) == 0];
[soldlv7, soldlh7] = solve(eqns, vars);

% disp(double(soldlv7));
% disp(double(soldlh7));
% disp('--------------------------');

syms lv7dd lh7dd;

vars = [lv7dd, lh7dd];

eqns = [l2*thdd2*sin(th2) - l2*(thd2^2)*cos(th2) + lv7dd*cos(thv7) + lh7dd*cos(thh7) == 0,...
        -l2*thdd2*cos(th2) + l2*(thd2^2)*sin(th2) + lv7dd*sin(thv7) + lh7dd*sin(thh7) == 0];
[solddlv7, solddlh7] = solve(eqns, vars);

% disp(double(solddlv7));
% disp(double(solddlh7));
% disp('--------------------------');

%% Vector Loop #3
syms th5 th8;

vars = [th5, th8];

eqns = [sollv7*cos(thv7) + sollh7*cos(thh7) + l3p*cos(sol3n+beta) + l5*cos(th5) - l8*cos(th8) == 0,...
        sollv7*sin(thv7) + sollh7*sin(thh7) + l3p*sin(sol3n+beta) + l5*sin(th5) - l8*sin(th8) == 0];
[sol5, sol8] = solve(eqns, vars);

% disp(double(sol5 * 180 / pi));
% disp(double(sol8 * 180 / pi));
% disp('--------------------------');

for i = 1:2
    if double(sol5(i) * 180 / pi) >= 90
        sol5n = double(sol5(i));
        sol8n = double(sol8(i));
    end
end

% disp(double(sol5n * 180 / pi));
% disp(double(sol8n * 180 / pi));
% disp(double(sol5n));
% disp(double(sol8n));
% disp('--------------------------');

syms thd5 thd8;

vars = [thd5, thd8];

eqns = [soldlv7*cos(thv7) + soldlh7*cos(thh7) - l3p*sold3*sin(sol3n+beta) - l5*thd5*sin(sol5n) + l8*thd8*sin(sol8n) == 0,...
        soldlv7*sin(thv7) + soldlh7*sin(thh7) + l3p*sold3*cos(sol3n+beta) + l5*thd5*cos(sol5n) - l8*thd8*cos(sol8n)];
[sold5, sold8] = solve(eqns, vars);

% disp(double(sold5));
% disp(double(sold8));
% disp('--------------------------');

syms thdd5 thdd8;

vars = [thdd5, thdd8];

eqns = [solddlv7*cos(thv7) + solddlh7*cos(thh7) - l3p*soldd3*sin(sol3n+beta) - l3p*(sold3^2)*cos(sol3n+beta) - l5*thdd5*sin(sol5n) - l5*(sold5^2)*cos(sol5n) + l8*thdd8*sin(sol8n) + l8*(sold8^2)*cos(sol8n) == 0,...
        solddlv7*sin(thv7) + solddlh7*sin(thh7) + l3p*soldd3*cos(sol3n+beta) - l3p*(sold3^2)*sin(sol3n+beta) + l5*thdd5*cos(sol5n) - l5*(sold5^2)*sin(sol5n) - l8*thdd8*cos(sol8n) + l8*(sold8^2)*sin(sol8n) == 0];
[soldd5, soldd8] = solve(eqns, vars);

% disp(double(soldd5));
% disp(double(soldd8));
% disp('--------------------------');

%% Final

th6 = sol8n + gamma;
thd6 = sold8;
thdd6 = soldd8;

% disp(double(th6 * 180 / pi));
% disp(double(th6));
% disp(double(thd6));
% disp(double(thdd6));
% disp('--------------------------');