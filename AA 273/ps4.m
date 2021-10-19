%% ps4  Problem Set 4
%
% Author: Tamas Kis
% Course: AA 273 - State Estimation and Filtering for Robotic Perception
% Last Update: 2021-08-23



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

% time parameters
dt = 0.001;     % time step [s]
tf = 10;        % simulation end time [s]
t = 0:dt:tf;    % time vector [s]
T = length(t);  % length of time vector

% principal inertia tensor [kg.m^2]
J = [1  0  0;
     0  5  0;
     0  0  5];

% threshold value for measurement model
c = 10;

% nonlinear dynamics and measurement equations
f = @(omega,tau,t) discrete_dynamics(omega,tau,t,J,dt);
g = @(omega,t) sat(omega,c);

% dynamics Jacobian and sensitivity matrix
A = @(omega,tau,t) dynamics_jacobian(omega,tau,t,J,dt);
C = @(omega,t) sensitivity_matrix(omega,t,c);

% state (n) and measurement (m) dimensions
n = 3;
m = 3;

% process (Q) and measurement (R) noise covariances
Q = 0.004*eye(n);
R = 0.01*eye(m);

% control input
tau = zeros(3,length(t));

% ground truth initial condition [rad/s]
omega0 = [10;0.1;0.1];

% initial condition structure for simulation
IC.x0_true = omega0;

% ground truth simulation
[omega,y] = simulate_nonlinear(f,g,Q,R,t,tau,IC,seed);



%% SIMULATION RESULTS

% ------------------
% x angular velocity
% ------------------

% initializes figure for wx
figure('position',pp.two_subplot_position);

% true state (wx)
subplot(1,2,1);
plot(t,omega(1,:),'color',pp.cardinal_red,'linewidth',pp.line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$\omega_{x}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize', pp.axis_font_size);
title('\textbf{True State}','interpreter','latex','fontsize',...
    pp.title_font_size)

% measurement (wx)
subplot(1,2,2);
plot(t(2:end),y(1,2:end),'color',pp.cardinal_red,'linewidth',pp.line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$\omega_{x}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize', pp.axis_font_size);
title('\textbf{Measurement}','interpreter','latex','fontsize',...
    pp.title_font_size)

% ------------------
% y angular velocity
% ------------------
% initializes figure for wy
figure('position',pp.two_subplot_position);

% true state (wy)
subplot(1,2,1);
plot(t,omega(2,:),'color',pp.cardinal_red,'linewidth',pp.line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$\omega_{y}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize', pp.axis_font_size);
title('\textbf{True State}','interpreter','latex','fontsize',...
    pp.title_font_size)

% measurement (wy)
subplot(1,2,2);
plot(t(2:end),y(2,2:end),'color',pp.cardinal_red,'linewidth',pp.line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$\omega_{y}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize', pp.axis_font_size);
title('\textbf{Measurement}','interpreter','latex','fontsize',...
    pp.title_font_size)

% ------------------
% z angular velocity
% ------------------

% initializes figure for wz
figure('position',pp.two_subplot_position);

% true state (wz)
subplot(1,2,1);
plot(t,omega(3,:),'color',pp.cardinal_red,'linewidth',pp.line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$\omega_{z}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize', pp.axis_font_size);
title('\textbf{True State}','interpreter','latex','fontsize',...
    pp.title_font_size)

% measurement (wz)
subplot(1,2,2);
plot(t(2:end),y(3,2:end),'color',pp.cardinal_red,'linewidth',pp.line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$\omega_{z}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize', pp.axis_font_size);
title('\textbf{Measurement}','interpreter','latex','fontsize',...
    pp.title_font_size)



%% FILTERING

% initial prior distribution (i.e. initial state estimate and covariance)
mu0 = [10;0;0];
Sigma0 = eye(3);

% runs extended Kalman filter
[mu,Sigma] = EKF(f,g,A,C,Q,R,t,tau,y,mu0,Sigma0);



%% FILTERING RESULTS

% initializes figure
figure('position',pp.three_subplot_position);

% true wx and its estimate
subplot(1,3,1);
hold on;
plot(t,omega(1,:),'linewidth',pp.line_width);
plot(t(2:end),mu(1,2:end),'linewidth',pp.line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$\omega_{x}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize',pp.axis_font_size);
legend('true state','estimated state','interpreter','latex','fontsize',...
    pp.legend_font_size);

% true wy and its estimate
subplot(1,3,2);
hold on;
plot(t,omega(2,:),'linewidth',pp.line_width);
plot(t(2:end),mu(2,2:end),'linewidth',pp.line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$\omega_{y}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize',pp.axis_font_size);
legend('true state','estimated state','interpreter','latex','fontsize',...
    pp.legend_font_size);

% true wz and its estimate
subplot(1,3,3);
hold on;
plot(t,omega(3,:),'linewidth',pp.line_width);
plot(t(2:end),mu(3,2:end),'linewidth',pp.line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$\omega_{z}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize',pp.axis_font_size);
legend('true state','estimated state','interpreter','latex','fontsize',...
    pp.legend_font_size);



%% ADDITIONAL FUNCTIONS

%==========================================================================
% discrete_dynamics  Discrete dynamics equation.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   omega       - (3×1 double) angular velocity (state vector) at current
%                 sample time [rad/s]
%   tau         - (3×1 double) control torque (control input) at current 
%                 sample time [N.m]
%   t           - (1×1 double) current time [s]
%   J           - (3×3 double) principal inertia tensor [kg.m^2]
%   dt          - (1×1 double) time step [s]
%
% -------
% OUTPUT:
% -------
%   omega_next  - (3×1 double) angular velocity (state vector) at next
%                 sample time [rad/s]
%
%==========================================================================
function omega_next = discrete_dynamics(omega,tau,t,J,dt)

    % unpacks angular vector
    wx = omega(1);
    wy = omega(2);
    wz = omega(3);
    
    % unpacks principal inertia tensor
    Jx = J(1,1);
    Jy = J(2,2);
    Jz = J(3,3);
    
    % evalutes f(x,u,t)
    omega_next = [wx+((Jy-Jz)*wy*wz*dt)/Jx;
                  wy+((Jz-Jx)*wz*wx*dt)/Jy;
                  wz+((Jx-Jy)*wx*wy*dt)/Jz]+[dt/Jx,dt/Jy,dt/Jz]*tau;
    
end



%==========================================================================
% satx  Element-wise saturation function (discrete measurement equation).
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (n×1 double) independent variable
%   c       - (1×1 double) threshold value for saturation
%
% -------
% OUTPUT:
% -------
%   satx    - (n×1 double) evaluation of sat(x)
%
% -----
% NOTE:
% -----
%   --> n = dimension of independent variable
%
%==========================================================================
function satx = sat(x,c)
    satx = zeros(size(x));
    for i = 1:length(x)
        if abs(x(i)) < c
            satx(i) = x(i);
        else
            satx(i) = c*sign(x(i));
        end
    end
end



%==========================================================================
% dynamics_Jacobian  Dynamics Jacobian.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   omega   - (3×1 double) angular velocity (state vector) [rad/s]
%   tau     - (3×1 double) control torque (control input) at current sample
%             time [N.m]
%   t       - (1×1 double) current time [s]
%   J       - (3×3 double) principal inertia tensor [kg.m^2]
%   dt      - (1×1 double) time step [s]
%
% -------
% OUTPUT:
% -------
%   A       - (3×3 double) dynamics Jacobian
%
%==========================================================================
function A = dynamics_jacobian(omega,tau,t,J,dt)

    % unpacks angular velocity vector [rad/s]
    wx = omega(1);
    wy = omega(2);
    wz = omega(3);
    
    % unpacks principal inertia tensor [kg.m^2]
    Jx = J(1,1);
    Jy = J(2,2);
    Jz = J(3,3);
    
    % assembles dynamics Jacobian
    A = [1                  (Jy-Jz)*wz*dt/Jx   (Jy-Jz)*wy*dt/Jx;
         (Jz-Jx)*wz*dt/Jy   1                  (Jz-Jx)*wx*dt/Jy;
         (Jx-Jy)*wy*dt/Jz   (Jx-Jy)*wx*dt/Jz   1];
    
end



%==========================================================================
% sensitivity_matrix  Sensitivity matrix.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   omega   - (3×1 double) angular velocity (state vector) [rad/s]
%   t       - (1×1 double) current time [s]
%   c       - (1×1 double) threshold value for measurement saturation 
%             [rad/s]
%
% -------
% OUTPUT:
% -------
%   C   	- (3×3 double) sensitivity matrix
%
%==========================================================================
function C = sensitivity_matrix(omega,t,c)

    % unpacks angular velocity
    wx = omega(1);
    wy = omega(2);
    wz = omega(3);
    
    % assembles sensitivity matrix
    C = [h(wx,c)   0         0;
         0         h(wy,c)   0;
         0         0         h(wz,c)];
    
end



%==========================================================================
% h  Derivative of the saturation function.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (1×1 double) independent variable
%   c       - (1×1 double) threshold value for saturation
%
% -------
% OUTPUT:
% -------
%   hx   	- (1×1 double) evaluation of h(x)
%
%==========================================================================
function hx = h(x,c)
    if abs(x) < c
        hx = 1;
    else
        hx = 0;
    end
end