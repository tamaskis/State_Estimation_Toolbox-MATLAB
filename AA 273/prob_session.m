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
PLOT_PARAMETERS



%% SIMULATION OF NONLINEAR DISCRETE TIME SYSTEM

% parameters
m = 1;
g0 = 10;
L = 1;

% initial condition
x0 = [pi/4;1];

% threshold value for measurement model
c = 10;

% dimension parameters
n = 2;
k = 1;

% noise covariances
Q = 0.0005*eye(n);
R = 0.4;

% time parameters
dt = 0.01; % time step [s]
tf = 10; % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% simulated noise
w = gaussian_random_sample(zeros(n,1),Q,T-1);
v = gaussian_random_sample(zeros(k,1),R,T-1);
w = zeros(size(w));

% control signal
u = zeros(1,length(t));

% array preallocation
x = zeros(n,T);
y = zeros(k,T);

% stores initial condition
x(:,1) = x0;

% functions for nonlinear dynamics and nonlinear measurements
f = @(x,u) x+dt*[x(2);-g0/L*sin(x(1))]+dt*[0;1/(m*L^2)]*u;
g = @(x) L*sin(x(1));

% simulation
for tt = 1:(T-1)
    
    % propagates state vector with noise
    x(:,tt+1) = f(x(:,tt),u(:,tt))+w(:,tt);
    
    % measurement with noise
    y(:,tt+1) = g(x(:,tt))+v(:,tt);
    
end



%% FILTERING

% initial state estimate
mu0 = [pi/5;1.1];
Sigma0 = 0.5*eye(2);

% dynamics and measurement Jacobians
A = @(x,u) [1,dt;-g0*(dt/L)*cos(x(1)),1];
C = @(x) [L*cos(x(1)),0];

% runs extended Kalman filter
[mu,Sigma] = EKF(f,g,A,C,Q,R,u,y,mu0,Sigma0);

% angle
figure('position',plot_position);
hold on;
plot(t,x(1,:),'linewidth',line_width);
plot(t(2:end),mu(1,2:end),'linewidth',line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\theta\;[\mathrm{rad}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true angle','estimated angle','interpreter','latex','fontsize',...
    legend_font_size);

% angular velocity
figure('position',plot_position);
hold on;
plot(t,x(2,:),'linewidth',line_width);
plot(t(2:end),mu(2,2:end),'linewidth',line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\dot{\theta}\;[\mathrm{rad/s}]$','interpreter','latex',...
    'fontsize',axis_font_size);
legend('true angular velocity','estimated angular velocity',...
    'interpreter','latex','fontsize',legend_font_size);