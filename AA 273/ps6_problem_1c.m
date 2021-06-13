%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Problem Set 6 Problem 1c

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

% seeds random number generators
rng(1);



%% SIMULATION OF NONLINEAR DISCRETE TIME SYSTEM

% features
m1 = [0;0];
m2 = [10;0];
m3 = [10;10];
m4 = [0;10];

% time parameters
dt = 0.1; % time step [s]
tf = 20; % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% dimension parameters
n = 3;
k = 4;

% noise covariances
Q = 0.1*eye(3)*dt;
R = 0.1*eye(4);

% initial state estimate (prior distribution)
mu0 = [0;0;0];
Sigma0 = 0.01*eye(3);

% simulated noise
w = gaussian_random_sample(zeros(n,1),Q,T-1);
v = gaussian_random_sample(zeros(k,1),R,T-1);

% control signal
u = [ones(size(t));sin(t)];

% array preallocation
x = zeros(n,T);
y = zeros(k,T);

% samples initial condition from the prior distribution
x(:,1) = gaussian_random_sample(mu0,Sigma0);

% functions for nonlinear dynamics and nonlinear measurements
f = @(x,u,t) f_discrete(x,u,t,dt);
g = @(x,t) g_discrete(x,t,m1,m2,m3,m4);

% simulation
for tt = 1:(T-1)
    
    % propagates state vector with noise
    x(:,tt+1) = f(x(:,tt),u(:,tt),tt)+w(:,tt);
    
    % measurement with noise
    y(:,tt+1) = g(x(:,tt),tt)+v(:,tt);
    
end

% dynamics and measurement Jacobians
A = @(x,u,t) dynamics_jacobian(x,u,t,dt);
C = @(x,t) measurement_jacobian(x,t,m1,m2,m3,m4);

% runs EKF
[mu,Sigma] = EKF_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,1:k);

% +/- 2-sigma bounds
[lower_bound,upper_bound] = sigma_bounds(mu,Sigma,2);

% initializes figure
figure('position',three_subplot_position);

% true x-position and its estimate
subplot(1,3,1);
hold on;
patch([t,fliplr(t)],[upper_bound(1,:),fliplr(lower_bound(1,:))],...
    matlab_light_red,'edgecolor','none','handlevisibility','off');
plot(t,x(1,:),'linewidth',line_width);
plot(t(2:end),mu(1,2:end),'linewidth',line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$p_{x}$','interpreter','latex','fontsize',axis_font_size);
legend('true $x$-position',...
    'estimated $x$-position (with $\pm2\sigma$ bound)','interpreter',...
    'latex','fontsize',legend_font_size);

% true y-position and its estimate
subplot(1,3,2);
hold on;
patch([t,fliplr(t)],[upper_bound(2,:),fliplr(lower_bound(2,:))],...
    matlab_light_red,'edgecolor','none','handlevisibility','off');
plot(t,x(2,:),'linewidth',line_width);
plot(t(2:end),mu(2,2:end),'linewidth',line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$p_{y}$','interpreter','latex','fontsize',axis_font_size);
legend('true $y$-position',...
    'estimated $y$-position (with $\pm2\sigma$ bound)','interpreter',...
    'latex','fontsize',legend_font_size);

% true heading angle and its estimate
subplot(1,3,3);
hold on;
patch([t,fliplr(t)],[upper_bound(3,:),fliplr(lower_bound(3,:))],...
    matlab_light_red,'edgecolor','none','handlevisibility','off');
plot(t,x(3,:),'linewidth',line_width);
plot(t(2:end),mu(3,2:end),'linewidth',line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\theta\;[\mathrm{rad}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true heading angle',...
    'estimated heading angle (with $\pm2\sigma$ bound)','interpreter',...
    'latex','fontsize',legend_font_size);



%% ADDITIONAL FUNCTIONS

%=========================================================================%
% Discrete-time nonlinear dynamics.
%=========================================================================%
%
% INPUTS:
%   x       state vector at current time step (3 x 1)
%   u       control input at current time step (2 x 1)
%   t       current iteration (corresponding to discrete time)
%   dt      time step [s]
%
% OUTPUTS:
%   f_eval  evaluation of f (produces state vector at next time step)
%
%=========================================================================%
function f_eval = f_discrete(x,u,t,dt)

    % unpacks state vector
    px = x(1);
    py = x(2);
    theta = x(3);
    
    % unpacks control vector
    nu = u(1);
    phi = u(2);
    
    % evalutes f(xt,ut)
    f_eval = [px+nu*cos(theta)*dt;py+nu*sin(theta)*dt;theta+phi*dt];
    
end



%=========================================================================%
% Discrete-time nonlinear measurement.
%=========================================================================%
%
% INPUTS:
%   x       state vector at current time step (3 x 1)
%   t       current iteration (corresponding to discrete time)
%   m1      feature 1 location (2 x 1)
%   m2      feature 2 location (2 x 1)
%   m3      feature 3 location (2 x 1)
%   m4      feature 4 location (2 x 1)
%
% OUTPUTS:
%   g_eval  evaluation of g (produces measurement at current time step)
%
%=========================================================================%
function g_eval = g_discrete(x,t,m1,m2,m3,m4)

    % unpacks state vector
    px = x(1);
    py = x(2);
    theta = x(3);
    
    % x and y components of features
    m1x = m1(1);
    m1y = m1(2);
    m2x = m2(1);
    m2y = m2(2);
    m3x = m3(1);
    m3y = m3(2);
    m4x = m4(1);
    m4y = m4(2);
    
    % evalutes g(xt)
    g_eval = [atan2(m1y-py,m1x-px)-theta;
              atan2(m2y-py,m2x-px)-theta;
              atan2(m3y-py,m3x-px)-theta;
              atan2(m4y-py,m4x-px)-theta];
                      
end



%=========================================================================%
% Dynamics Jacobian.
%=========================================================================%
%
% INPUTS:
%   x       state vector (3 x 1)
%   u       control input (2 x 1)
%   t       current iteration (corresponding to discrete time)
%   dt      time step [s]
%
% OUTPUTS:
%   A       dynamics Jacobian
%
%=========================================================================%
function A = dynamics_jacobian(x,u,t,dt)

    % extracts "theta" from state vector
    theta = x(3);
    
    % extracts "nu" from control input
    nu = u(1);
    
    % assembles dynamics Jacobian
    A = [1   0   -nu*sin(theta)*dt;
         0   1    nu*cos(theta)*dt;
         0   0    1];
    
end



%=========================================================================%
% Measurement Jacobian.
%=========================================================================%
%
% INPUTS:
%   x       state vector (3 x 1)
%   t       current iteration (corresponding to discrete time)
%   m1      feature 1 location (2 x 1)
%   m2      feature 2 location (2 x 1)
%   m3      feature 3 location (2 x 1)
%   m4      feature 4 location (2 x 1)
%
% OUTPUTS:
%   C       measurement Jacobian
%
%=========================================================================%
function C = measurement_jacobian(x,t,m1,m2,m3,m4)
    
    % extracts "px" and "py" from state vector
    px = x(1);
    py = x(2);
    
    % x and y components of features
    m1x = m1(1);
    m1y = m1(2);
    m2x = m2(1);
    m2y = m2(2);
    m3x = m3(1);
    m3y = m3(2);
    m4x = m4(1);
    m4y = m4(2);
    
    % assembles measurement Jacobian
    C = [(m1y-py)/((m1x-px)^2+(m1y-py)^2)   (px-m1x)/((m1x-px)^2+(m1y-py)^2)   -1;
         (m2y-py)/((m2x-px)^2+(m2y-py)^2)   (px-m2x)/((m2x-px)^2+(m2y-py)^2)   -1;
         (m3y-py)/((m3x-px)^2+(m3y-py)^2)   (px-m3x)/((m3x-px)^2+(m3y-py)^2)   -1;
         (m4y-py)/((m4x-px)^2+(m4y-py)^2)   (px-m4x)/((m4x-px)^2+(m4y-py)^2)   -1];

end