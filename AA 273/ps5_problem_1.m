%% ps5_problem_1
% Problem Set 5, Problem 1
% AA 273 - State Estimation and Filtering for Robotic Perception
%
% Author: Tamas Kis
% Last Update: 2021-08-18



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;

% adds path to root directory and all subdirectories
addpath(genpath("../"));

% loads plot parameters
pp = PLOT_PARAMETERS;

% random seed
seed = 2;



%% SIMULATION

% roll mode parameters
a = 0.5;
b = 2;

% dimension parameters
n = 3;
k = 1;

% noise covariances
Q = [0.1,0,0;0,0,0;0,0,0];
R = 0.1;

% time parameters
dt = 0.1;   % time step [s]
tf = 20;    % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% control input
u = sin(t);

% samples initial condition for roll rate [rad/s]
p0 = gaussian_random_sample(0,0.1);

% initial state
x0 = [p0;a;b];

% functions for nonlinear dynamics and nonlinear measurements
f = @(x,u,t) discrete_dynamics(x,u,t,dt);
g = @(x,t) discrete_measurement(x,t);

% initial condition structure for simulation
IC.x0_true = x0;

% ground truth simulation
[x,y] = simulate_nonlinear(f,g,Q,R,t,u,IC,seed);



%% SIMULATION RESULTS

% initializes figure for roll rate
figure('position',pp.two_subplot_position);

% true roll rate
subplot(1,2,1);
plot(t,x(1,:),'color',pp.cardinal_red,'linewidth',pp.line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$p\;[\mathrm{rad/s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
title('\textbf{True Roll Rate}','interpreter','latex','fontsize',...
    pp.title_font_size)

% measured roll rate
subplot(1,2,2);
plot(t(2:end),y(1,2:end),'color',pp.cardinal_red,'linewidth',...
    pp.line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$p\;[\mathrm{rad/s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
title('\textbf{Measured Roll Rate}','interpreter','latex','fontsize',...
    pp.title_font_size)



%% FILTERING

% initial state estimate
mu0 = [0;0;0];
Sigma0 = 10*eye(3);

% dynamics and measurement Jacobians
A = @(x,u,t) dynamics_jacobian(x,u,t,dt);
C = @(x,t) measurement_jacobian(x,t);

% redefines Q for EKF
Q = [0.1,0,0;0,0.01,0;0,0,0.01];

% runs EKF
[mu,Sigma] = EKF(f,g,A,C,Q,R,t,u,y,mu0,Sigma0);

% 2-sigma bounds
[mu_minus,mu_plus] = sigma_bounds(mu,Sigma,2);



%% FILTER RESULTS

% initializes figure
figure('position',pp.three_subplot_position);

% true roll rate and its estimate
subplot(1,3,1);
hold on;
patch([t,fliplr(t)],[mu_plus(1,:),fliplr(mu_minus(1,:))],...
    pp.matlab_light_red,'edgecolor','none','handlevisibility','off');
plot(t,x(1,:),'linewidth',pp.line_width);
plot(t(2:end),mu(1,2:end),'linewidth',pp.line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$p\;[\mathrm{rad/s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
legend('true roll rate',...
    'estimated roll rate (with $\pm2\sigma$ bounds)','interpreter',...
    'latex','fontsize',pp.legend_font_size);

% true a and its estimate
subplot(1,3,2);
hold on;
patch([t,fliplr(t)],[mu_plus(2,:),fliplr(mu_minus(2,:))],...
    pp.matlab_light_red,'edgecolor','none','handlevisibility','off');
plot(t,x(2,:),'linewidth',pp.line_width);
plot(t(2:end),mu(2,2:end),'linewidth',pp.line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$a$','interpreter','latex','fontsize',pp.axis_font_size);
legend('true $a$','estimated $a$ (with $\pm2\sigma$ bounds)',...
    'interpreter','latex','fontsize',pp.legend_font_size);

% true b and its estimate
subplot(1,3,3);
hold on;
patch([t,fliplr(t)],[mu_plus(3,:),fliplr(mu_minus(3,:))],...
    pp.matlab_light_red,'edgecolor','none','handlevisibility','off');
plot(t,x(3,:),'linewidth',pp.line_width);
plot(t(2:end),mu(3,2:end),'linewidth',pp.line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$b$','interpreter','latex','fontsize',pp.axis_font_size);
legend('true $b$','estimated $b$ (with $\pm2\sigma$ bounds)',...
    'interpreter','latex','fontsize',pp.legend_font_size);



%% ADDITIONAL FUNCTIONS

%==========================================================================
% discrete_dynamics  Discrete dynamics equation.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (3×1 double) state vector at current sample time
%   u       - (1×1 double) control input at current sample time
%   t       - (1×1 double) current time [s]
%   dt      - (1×1 double) time step [s]
%
% -------
% OUTPUT:
% -------
%   x_next  - (3×1 double) state vector at next sample time
%
%==========================================================================
function x_next = discrete_dynamics(x,u,t,dt)

    % unpacks state vector
    p = x(1);
    a = x(2);
    b = x(3);
    
    % evalutes f(x,u)
    x_next = [p+(-a*p+b*u)*dt;a;b];
    
end



%==========================================================================
% discrete_measurement  Discrete measurement equation.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (3×1 double) state vector
%   t       - (1×1 double) time [s]
%
% -------
% OUTPUT:
% -------
%   y       - (1×1 double) measurement
%
%==========================================================================
function y = discrete_measurement(x,t)

    % extracts "p" from state vector
    p = x(1);
    
    % evalutes g(x)
    y = p;
    
end



%==========================================================================
% dynamics_jacobian  Dynamics Jacobian.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x	- (3×1 double) state vector
%   u   - (1×1 double) control input
%   t 	- (1×1 double) time [s]
%   dt  - (1×1 double) time step [s]
%
% -------
% OUTPUT:
% -------
%   A   - (3×3 double) dynamics Jacobian
%
%==========================================================================
function A = dynamics_jacobian(x,u,t,dt)

    % extacts "p" and "a" from state vector
    p = x(1);
    a = x(2);
    
    % assembles dynamics Jacobian
    A = [1-a*dt   -p*dt   u*dt;
         0         1      0;
         0         0      1];
    
end



%==========================================================================
% measurement_jacobian  Measurement Jacobian.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x   - (3×1 double) state vector
%   t 	- (1×1 double) time [s]
%
% -------
% OUTPUT:
% -------
%   C   - (3×1 double) measurement Jacobian
%
%==========================================================================
function C = measurement_jacobian(x,t)
    C = [1,0,0];
end