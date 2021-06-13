%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Problem Set 4 Problem 3

%-------------------------------------------------------------------------%



%% SCRIPT SETUP

% clears variables and command window, closes all figures
clear;
clc;
close all;

% adds path to Estimation Toolbox
addpath("..");

% loads plot parameters
PLOT_PARAMETERS;



%% SIMULATION OF NONLINEAR DISCRETE TIME SYSTEM

% principal inertia tensor
J = [1,0,0;0,5,0;0,0,5];

% initial condition
omega0 = [10;0.1;0.1]; % angular velocity [rad/s]

% threshold value for measurement model
c = 10;

% dimension parameters
n = 3;
k = 3;

% noise covariances
Q = 0.004*eye(n);
R = 0.01*eye(k);

% time parameters
dt = 0.001; % time step [s]
tf = 10; % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% simulated noise
w = gaussian_random_sample(zeros(n,1),Q,T-1);
v = gaussian_random_sample(zeros(k,1),R,T-1);

% control signal
tau = zeros(3,length(t));

% array preallocation
omega = zeros(n,T);
y = zeros(k,T);

% stores initial condition
omega(:,1) = omega0;

% functions for nonlinear dynamics and nonlinear measurements
f = @(omega,tau) f_discrete(omega,tau,J,dt);
g = @(omega) sat(omega,c);

% simulation
for tt = 1:(T-1)
    
    % propagates state vector with noise
    omega(:,tt+1) = f(omega(:,tt),tau(:,tt))+w(:,tt);
    
    % measurement with noise
    y(:,tt+1) = g(omega(:,tt))+v(:,tt);
    
end

% initializes figure for wx
figure('position',two_subplot_position);

% true state (wx)
subplot(1,2,1);
plot(t,omega(1,:),'color',cardinal_red,'linewidth',line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\omega_{x}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize', axis_font_size);
title('\textbf{True State}','interpreter','latex','fontsize',...
    title_font_size)

% measurement (wx)
subplot(1,2,2);
plot(t(2:end),y(1,2:end),'color',cardinal_red,'linewidth',line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\omega_{x}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize', axis_font_size);
title('\textbf{Measurement}','interpreter','latex','fontsize',...
    title_font_size)

% initializes figure for wy
figure('position',two_subplot_position);

% true state (wy)
subplot(1,2,1);
plot(t,omega(2,:),'color',cardinal_red,'linewidth',line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\omega_{y}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize', axis_font_size);
title('\textbf{True State}','interpreter','latex','fontsize',...
    title_font_size)

% measurement (wy)
subplot(1,2,2);
plot(t(2:end),y(2,2:end),'color',cardinal_red,'linewidth',line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\omega_{y}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize', axis_font_size);
title('\textbf{Measurement}','interpreter','latex','fontsize',...
    title_font_size)

% initializes figure for wz
figure('position',two_subplot_position);

% true state (wz)
subplot(1,2,1);
plot(t,omega(3,:),'color',cardinal_red,'linewidth',line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\omega_{z}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize', axis_font_size);
title('\textbf{True State}','interpreter','latex','fontsize',...
    title_font_size)

% measurement (wz)
subplot(1,2,2);
plot(t(2:end),y(3,2:end),'color',cardinal_red,'linewidth',line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\omega_{z}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize', axis_font_size);
title('\textbf{Measurement}','interpreter','latex','fontsize',...
    title_font_size)



%% FILTERING

% initial state estimate
mu0 = [10;0;0];
Sigma0 = eye(3);

% dynamics and measurement Jacobians
A = @(omega,tau) dynamics_jacobian(omega,tau,J,dt);
C = @(omega) measurement_jacobian(omega,c);

% runs extended Kalman filter
[mu,Sigma] = EKF(f,g,A,C,Q,R,tau,y,mu0,Sigma0);

% initializes figure
figure('position',three_subplot_position);

% true wx and its estimate
subplot(1,3,1);
hold on;
plot(t,omega(1,:),'linewidth',line_width);
plot(t(2:end),mu(1,2:end),'linewidth',line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\omega_{x}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize',axis_font_size);
legend('true state','estimated state','interpreter','latex','fontsize',...
    legend_font_size);

% true wy and its estimate
subplot(1,3,2);
hold on;
plot(t,omega(2,:),'linewidth',line_width);
plot(t(2:end),mu(2,2:end),'linewidth',line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\omega_{y}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize',axis_font_size);
legend('true state','estimated state','interpreter','latex','fontsize',...
    legend_font_size);

% true wz and its estimate
subplot(1,3,3);
hold on;
plot(t,omega(3,:),'linewidth',line_width);
plot(t(2:end),mu(3,2:end),'linewidth',line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\omega_{z}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize',axis_font_size);
legend('true state','estimated state','interpreter','latex','fontsize',...
    legend_font_size);



%% ADDITIONAL FUNCTIONS

%=========================================================================%
% Discrete-time nonlinear dynamics.
%=========================================================================%
% INPUT: omega - angular velocity (state vector) [rad/s]
%        tau - torque (control input) [N.m]
%        J - principal inertia tensor [kg.m^2]
%        dt - time step [s]
% OUTPUT: omega_t1 - angular velocity at next time step [rad/s]
function omega_t1 = f_discrete(omega_t,tau_t,J,dt)

    % unpacks angular vector
    wx = omega_t(1);
    wy = omega_t(2);
    wz = omega_t(3);
    
    % unpacks principal inertia tensor
    Jx = J(1,1);
    Jy = J(2,2);
    Jz = J(3,3);
    
    % evalutes f(omega_t,tau_t)
    omega_t1 = [wx+((Jy-Jz)*wy*wz*dt)/Jx;wy+((Jz-Jx)*wz*wx*dt)/Jy;...
        wz+((Jx-Jy)*wx*wy*dt)/Jz]+[dt/Jx,dt/Jy,dt/Jz]*tau_t;
    
end



%=========================================================================%
% Dynamics Jacobian.
%=========================================================================%
% INPUT: omega - angular velocity (state vector) [rad/s]
%        tau - torque (control input) [N.m]
%        J - principal inertia tensor [kg.m^2]
%        dt - time step [s]
% OUTPUT: A - dynamics Jacobian
function A = dynamics_jacobian(omega,tau,J,dt)

    % unpacks angular vector
    wx = omega(1);
    wy = omega(2);
    wz = omega(3);
    
    % unpacks principal inertia tensor
    Jx = J(1,1);
    Jy = J(2,2);
    Jz = J(3,3);
    
    % assembles dynamics Jacobian
    A = [1                  (Jy-Jz)*wz*dt/Jx   (Jy-Jz)*wy*dt/Jx;
         (Jz-Jx)*wz*dt/Jy   1                  (Jz-Jx)*wx*dt/Jy;
         (Jx-Jy)*wy*dt/Jz   (Jx-Jy)*wx*dt/Jz   1];
    
end



%=========================================================================%
% Measurement Jacobian.
%=========================================================================%
% INPUT: omega - angular velocity (state vector) [rad/s]
%        h - function handle for derivative of saturation function
%        c - threshold value for measurement saturation [rad/s]
% OUTPUT: C - measurement Jacobian
function C = measurement_jacobian(omega,c)

    % unpacks angular vector
    wx = omega(1);
    wy = omega(2);
    wz = omega(3);
    
    % assembles measurement Jacobian
    C = [h(wx,c)   0         0;
         0         h(wy,c)   0;
         0         0         h(wz,c)];
    
end



%=========================================================================%
% Element-wise saturation function.
%=========================================================================%
% INPUT: x - independent variable
%        c - threshold value for saturation
% OUTPUT: satx - evaluation of sat(x)
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



%=========================================================================%
% Derivative of saturation function.
%=========================================================================%
% INPUT: x - independent variable
%        c - threshold value
% OUTPUT: hx - evaluation of h(x)
function hx = h(x,c)
    if abs(x) < c
        hx = 1;
    else
        hx = 0;
    end
end