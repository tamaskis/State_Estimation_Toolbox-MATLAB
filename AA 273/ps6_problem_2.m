%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Problem Set 6 Problem 2

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
rng(2);



%% SIMULATION OF NONLINEAR DISCRETE TIME SYSTEM

% time parameters
dt = 0.1; % time step [s]
tf = 50; % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% dimension parameters
n = 8;
k = 4;

% noise covariances
Q = 0.01*eye(n)*dt;
R = 0.1*eye(k);

% true feature locations
m1 = [0;0];
m2 = [10;0];
m3 = [10;10];
m4 = [0;10];

% initial state estimate (prior distribution)
mu0 = [m1;m2;m3;m4]+5;
Sigma0 = 0.01*eye(length(mu0));

% simulated noise
w = gaussian_random_sample(zeros(n,1),Q,T-1);
v = gaussian_random_sample(zeros(k,1),R,T-1);

% control signal
u = [ones(size(t));sin(t)];

% simulates robot trajectory
px = zeros(1,T);
py = zeros(1,T);
theta = zeros(1,T);
for tt = 1:(T-1)
    nu = u(1,tt);
    phi = u(2,tt);
    px(tt+1) = px(tt)+nu*cos(theta(tt))*dt+gaussian_random_sample(0,0.1*...
        dt);
    py(tt+1) = py(tt)+nu*sin(theta(tt))*dt+gaussian_random_sample(0,0.1*...
        dt);
    theta(tt+1) = theta(tt)+phi*dt+gaussian_random_sample(0,0.1*dt);
end

% array preallocation
x = zeros(n,T);
y = zeros(k,T);

% samples initial condition from the prior distribution
x(:,1) = gaussian_random_sample(mu0,Sigma0);

% functions for nonlinear dynamics and nonlinear measurements
f = @(x,u,t) f_discrete(x,u,t);
g = @(x,t) g_discrete(x,t,px,py,theta);

% simulation (using true feature locations)
for tt = 1:(T-1)
    
    % propagates state vector with noise
    x(:,tt+1) = f([m1;m2;m3;m4],u(:,tt),tt)+w(:,tt);
    
    % measurement with noise
    y(:,tt+1) = g([m1;m2;m3;m4],tt)+v(:,tt);
    
end



%% FILTERING

% dynamics and measurement Jacobians
A = @(x,u,t) dynamics_jacobian(x,u,dt);
C = @(x,t) measurement_jacobian(x,t,px,py);

% runs EKF
[mu,Sigma] = EKF_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,1:k);

% extracts time history of feature locations
m1_ekf = mu(1:2,2:end);
m2_ekf = mu(3:4,2:end);
m3_ekf = mu(5:6,2:end);
m4_ekf = mu(7:8,2:end);

% feature locations
figure('position',plot_position);
hold on;
plot(px,py,'linewidth',line_width,'color',cardinal_red);
scatter([m1(1),m2(1),m3(1),m4(1)],[m1(2),m2(2),m3(2),m4(2)],50,...
    'markeredgecolor','k','MarkerFaceColor','g','linewidth',1);
scatter([mu0(1),mu0(3),mu0(5),mu0(7)],[mu0(2),mu0(4),mu0(6),mu0(8)],50,...
    'markeredgecolor','k','MarkerFaceColor','r','linewidth',1);
scatter(m1_ekf(1,:),m1_ekf(2,:),5,matlab_blue,'filled');
scatter(m2_ekf(1,:),m2_ekf(2,:),5,matlab_red,'filled');
scatter(m3_ekf(1,:),m3_ekf(2,:),5,matlab_yellow,'filled');
scatter(m4_ekf(1,:),m4_ekf(2,:),5,matlab_purple,'filled');
hold off;
grid on;
axis equal;
xlim([-5,30]);
ylim([-5,30]);
xlabel('$p_{x}$','interpreter','latex','fontsize',axis_font_size);
ylabel('$p_{y}$','interpreter','latex','fontsize',axis_font_size);
legend('robot trajectory','true feature locations',...
    'initial guess of feature locations','feature 1 estimate',...
    'feature 2 estimate','feature 3 estimate','feature 4 estimate',...
    'interpreter','latex','fontsize',legend_font_size,'location',...
    'northeast');



%% ADDITIONAL FUNCTIONS

%=========================================================================%
% Discrete-time nonlinear dynamics.
%=========================================================================%
%
% INPUTS:
%   x       state vector at current time step (3 x 1)
%   u       control input at current time step (2 x 1)
%   t       current iteration (corresponding to discrete time)
%
% OUTPUTS:
%   f_eval  evaluation of f (produces state vector at next time step)
%
%=========================================================================%
function f_eval = f_discrete(x,u,t)
    f_eval = x;
end



%=========================================================================%
% Discrete-time nonlinear measurement.
%=========================================================================%
%
% INPUTS:
%   x       state vector at current time step (8 x 1)
%   t       current iteration (corresponding to discrete time)
%   px      x-position (1 x T)
%   py      y-position (1 x T)
%   theta   heading angle (1 x T) [rad]
%
% OUTPUTS:
%   g_eval  evaluation of g (produces measurement at current time step)
%
%=========================================================================%
function g_eval = g_discrete(x,t,px,py,theta)
    
    % unpacks state vector
    m1x = x(1);
    m1y = x(2);
    m2x = x(3);
    m2y = x(4);
    m3x = x(5);
    m3y = x(6);
    m4x = x(7);
    m4y = x(8);
    
    % overrides px, py, and theta with their values at time step t
    px = px(t);
    py = py(t);
    theta = theta(t);
    
    % evalutes g(xt,ut)
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
%   x       state vector (8 x 1)
%   u       control input (2 x 1)
%   t       current iteration (corresponding to discrete time)
%
% OUTPUTS:
%   A       dynamics Jacobian
%
%=========================================================================%
function A = dynamics_jacobian(x,u,t)
    A = eye(8);
end



%=========================================================================%
% Measurement Jacobian.
%=========================================================================%
%
% INPUTS:
%   x       state vector (8 x 1)
%   t       current iteration (corresponding to discrete time)
%   px      x-position (1 x T)
%   py      y-position (1 x T)
%
% OUTPUTS:
%   C       measurement Jacobian
%
%=========================================================================%
function C = measurement_jacobian(x,t,px,py)
    
    % creates "mx" and "my" vectors from state vector
    mx = x([1,3,5,7]);
    my = x([2,4,6,8]);
    
    % overrides px and py with their values at time step t
    px = px(t);
    py = py(t);
    
    % preallocates measurement Jacobian
    C = zeros(4,8);
    
    % assembles measurement Jacobian
    for i = 1:4
        for j = 1:4
            if i == j
                C(i,2*j-1) = (py-my(j))/((mx(j)-px)^2+(my(j)-py)^2);
                C(i,2*j) = (mx(j)-px)/((mx(j)-px)^2+(my(j)-py)^2);
            end
        end
    end
    
end