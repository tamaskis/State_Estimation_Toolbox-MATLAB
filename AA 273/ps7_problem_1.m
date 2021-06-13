%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Problem Set 7 Problem 1

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
dt = 1; % time step [s]
tf = 50; % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% dimension parameters
n = 6;
k = 6;

% noise covariances
Q = @(j) process_noise_covariance(j);
R = @(i) measurement_noise_covariance(i);

% initial state estimate (prior distribution)
mu0 = [0;0;25;0;0;25];
mu0 = [1;1;2;2;3;3];
Sigma0 = 0.1*eye(length(mu0));

% simulated noise
w = gaussian_random_sample(zeros(n,1),Q(1),T-1);
v = gaussian_random_sample(zeros(k,1),R(1),T-1);

% system, input, and measurement matrices
A = @(x,u,t,j) system_matrix(x,u,t,j);
B = @(x,u,t,j) input_matrix(x,u,t,j);
C = @(x,t,i) measurement_matrix(x,t,i);

% control vector
u = [cos(0.1*t);sin(0.1*t);-cos(0.2*t);sin(0.2*t);cos(0.1*t);sin(0.2*t)];

% array preallocation
x = zeros(n,T);
y = zeros(k,T);

% samples initial condition from initial prior distribution
x(:,1) = gaussian_random_sample(mu0,Sigma0,1);

% simulation
for tt = 1:(T-1)
    
    % propagates state vector with noise
    x(:,tt+1) = A(x(:,tt),u(:,tt),t,1)*x(:,tt)+B(x(:,tt),u(:,tt),t,1)*...
        u(:,tt)+w(:,tt);
    
    % measurement with noise
    y(:,tt+1) = C(x(:,tt),t,1)*x(:,tt+1)+v(:,tt);
    
end

% true target trajectories
x1 = x(1,:);
y1 = x(2,:);
x2 = x(3,:);
y2 = x(4,:);
x3 = x(5,:);
y3 = x(6,:);



%% FILTERING

% number of Gaussian components in dynamics model
M = 1;

% number of Gaussian components in sensor model
L = 3;

% number of Gaussian components to retain during pruning
Z = 3;

beta = [1/3;1/3;1/3];
gamma = ones(6,1)/6;

% runs MHKF
[mu,Sigma,alpha] = MHKF(A,B,C,Q,R,u,y,{mu0},{Sigma0},beta,gamma,M,L,Z);

% extracts results for various mixture components
mu1 = mu(:,2:end,1);
mu2 = mu(:,2:end,2);
mu3 = mu(:,2:end,3);
Sigma1 = Sigma(:,:,2:end,1);
Sigma2 = Sigma(:,:,2:end,2);
Sigma3 = Sigma(:,:,2:end,3);
alpha1 = alpha(1,:);
alpha2 = alpha(2,:);
alpha3 = alpha(3,:);

% mixture component 1
figure('position',plot_position);
hold on;
for i = 1:(length(t)-1)
    [x_ell,y_ell] = error_ellipse(mu1(1:2,i),Sigma1(1:2,1:2,i),0.95);
    patch(x_ell,y_ell,matlab_light_blue,'edgecolor','none',...
        'handlevisibility','off');
    [x_ell,y_ell] = error_ellipse(mu1(3:4,i),Sigma1(3:4,3:4,i),0.95);
    patch(x_ell,y_ell,matlab_light_red,'edgecolor','none',...
        'handlevisibility','off');
    [x_ell,y_ell] = error_ellipse(mu1(5:6,i),Sigma1(5:6,5:6,i),0.95);
    patch(x_ell,y_ell,matlab_light_yellow,'edgecolor','none',...
        'handlevisibility','off');
end
plot(mu1(1,:),mu1(2,:),'linewidth',line_width,'color',matlab_blue);
plot(x1,y1,'linewidth',line_width,'linestyle',':','color',matlab_blue);
plot(mu1(3,:),mu1(4,:),'linewidth',line_width,'color',matlab_red);
plot(x2,y2,'linewidth',line_width,'linestyle',':','color',matlab_red);
plot(mu1(5,:),mu1(6,:),'linewidth',line_width,'color',matlab_yellow);
plot(x3,y3,'linewidth',line_width,'linestyle',':','color',matlab_yellow);
hold off;
grid on;
xlabel('$x$','interpreter','latex','fontsize',axis_font_size);
ylabel('$y$','interpreter','latex','fontsize',axis_font_size);
title('\textbf{Mixture Component 1}','interpreter','latex','fontsize',...
    title_font_size);
legend('target 1 estimated trajectory','target 1 true trajectory',...
    'target 2 estimated trajectory','target 2 true trajectory',...
    'target 3 estimated trajectory','target 3 true trajectory',...
    'interpreter','latex','fontsize',legend_font_size,'location',...
    'northeast');

% mixture component 2
figure('position',plot_position);
hold on;
for i = 1:(length(t)-1)
    [x_ell,y_ell] = error_ellipse(mu2(1:2,i),Sigma2(1:2,1:2,i),0.95);
    patch(x_ell,y_ell,matlab_light_blue,'edgecolor','none',...
        'handlevisibility','off');
    [x_ell,y_ell] = error_ellipse(mu2(3:4,i),Sigma2(3:4,3:4,i),0.95);
    patch(x_ell,y_ell,matlab_light_red,'edgecolor','none',...
        'handlevisibility','off');
    [x_ell,y_ell] = error_ellipse(mu2(5:6,i),Sigma2(5:6,5:6,i),0.95);
    patch(x_ell,y_ell,matlab_light_yellow,'edgecolor','none',...
        'handlevisibility','off');
end
plot(mu2(1,:),mu2(2,:),'linewidth',line_width,'color',matlab_blue);
plot(x1,y1,'linewidth',line_width,'linestyle',':','color',matlab_blue);
plot(mu2(3,:),mu2(4,:),'linewidth',line_width,'color',matlab_red);
plot(x2,y2,'linewidth',line_width,'linestyle',':','color',matlab_red);
plot(mu2(5,:),mu2(6,:),'linewidth',line_width,'color',matlab_yellow);
plot(x3,y3,'linewidth',line_width,'linestyle',':','color',matlab_yellow);
hold off;
grid on;
xlabel('$x$','interpreter','latex','fontsize',axis_font_size);
ylabel('$y$','interpreter','latex','fontsize',axis_font_size);
title('\textbf{Mixture Component 2}','interpreter','latex','fontsize',...
    title_font_size);
legend('target 1 estimated trajectory','target 1 true trajectory',...
    'target 2 estimated trajectory','target 2 true trajectory',...
    'target 3 estimated trajectory','target 3 true trajectory',...
    'interpreter','latex','fontsize',legend_font_size,'location',...
    'northeast');

% mixture component 3
figure('position',plot_position);
hold on;
for i = 1:(length(t)-1)
    [x_ell,y_ell] = error_ellipse(mu3(1:2,i),Sigma3(1:2,1:2,i),0.95);
    patch(x_ell,y_ell,matlab_light_blue,'edgecolor','none',...
        'handlevisibility','off');
    [x_ell,y_ell] = error_ellipse(mu3(3:4,i),Sigma3(3:4,3:4,i),0.95);
    patch(x_ell,y_ell,matlab_light_red,'edgecolor','none',...
        'handlevisibility','off');
    [x_ell,y_ell] = error_ellipse(mu3(5:6,i),Sigma3(5:6,5:6,i),0.95);
    patch(x_ell,y_ell,matlab_light_yellow,'edgecolor','none',...
        'handlevisibility','off');
end
plot(mu3(1,:),mu3(2,:),'linewidth',line_width,'color',matlab_blue);
plot(x1,y1,'linewidth',line_width,'linestyle',':','color',matlab_blue);
plot(mu3(3,:),mu3(4,:),'linewidth',line_width,'color',matlab_red);
plot(x2,y2,'linewidth',line_width,'linestyle',':','color',matlab_red);
plot(mu3(5,:),mu3(6,:),'linewidth',line_width,'color',matlab_yellow);
plot(x3,y3,'linewidth',line_width,'linestyle',':','color',matlab_yellow);
hold off;
grid on;
xlabel('$x$','interpreter','latex','fontsize',axis_font_size);
ylabel('$y$','interpreter','latex','fontsize',axis_font_size);
title('\textbf{Mixture Component 3}','interpreter','latex','fontsize',...
    title_font_size);
legend('target 1 estimated trajectory','target 1 true trajectory',...
    'target 2 estimated trajectory','target 2 true trajectory',...
    'target 3 estimated trajectory','target 3 true trajectory',...
    'interpreter','latex','fontsize',legend_font_size,'location',...
    'northeast');

% mixture weights
figure('position',plot_position);
hold on;
plot(t,alpha1,'linewidth',line_width,'color',matlab_blue);
plot(t,alpha2,'linewidth',line_width,'color',matlab_red);
plot(t,alpha3,'linewidth',line_width,'color',matlab_yellow);
hold off;
grid on;
xlabel('$t$','interpreter','latex','fontsize',axis_font_size);
ylabel('$\alpha_{t|t}^{i}$','interpreter','latex','fontsize',...
    axis_font_size);
title('\textbf{Mixture Weights}','interpreter','latex','fontsize',...
    title_font_size);
legend('mixture component 1','mixture component 2',...
    'mixture component 3','interpreter','latex','fontsize',...
    legend_font_size,'location','northeast');



%% ADDITIONAL FUNCTIONS

%==========================================================================
% System (plant) matrix.
%==========================================================================
%
% INPUTS:
%   x       state vector (6 x 1)
%   u       control input (6 x 1)
%   t       current iteration (corresponding to discrete time)
%   j       index of dynamics equation
%
% OUTPUTS:
%   A       system (plant) matrix
%
%==========================================================================
function A = system_matrix(x,u,t,j)
    A = eye(6);
end



%==========================================================================
% Input (control) matrix.
%==========================================================================
%
% INPUTS:
%   x       state vector (6 x 1)
%   u       control input (6 x 1)
%   t       current iteration (corresponding to discrete time)
%   j       index of dynamics equation
%
% OUTPUTS:
%   B       input (control) matrix
%
%==========================================================================
function B = input_matrix(x,u,t,j)
    B = eye(6);
end



%==========================================================================
% Measurement matrix.
%==========================================================================
%
% INPUTS:
%   x       state vector (6 x 1)
%   t       current iteration (corresponding to discrete time)
%   i       index of measurement equation
%
% OUTPUTS:
%   C       measurement matrix (6 x 6)
%
%==========================================================================
function C = measurement_matrix(x,t,i)
    if i == 1
        C = [eye(2)     zeros(2)   zeros(2);
             zeros(2)   eye(2)     zeros(2);
             zeros(2)   zeros(2)   eye(2)];
    elseif i == 2
        C = [eye(2)     zeros(2)   zeros(2);
             zeros(2)   zeros(2)   eye(2);
             zeros(2)   eye(2)     zeros(2)];
    elseif i == 3
        C = [zeros(2)   eye(2)     zeros(2);
             eye(2)     zeros(2)   zeros(2);
             zeros(2)   zeros(2)   eye(2)];
    elseif i == 4
        C = [zeros(2)   eye(2)     zeros(2);
             zeros(2)   zeros(2)   eye(2);
             eye(2)     zeros(2)   zeros(2)];
    elseif i == 5
        C = [zeros(2)   zeros(2)   eye(2);
             eye(2)     zeros(2)   zeros(2);
             zeros(2)   eye(2)     zeros(2)];
    elseif i == 6
        C = [zeros(2)   zeros(2)   eye(2);
             zeros(2)   eye(2)     zeros(2);
             eye(2)     zeros(2)   zeros(2)];
    end
end



%==========================================================================
% Process noise covariance.
%==========================================================================
%
% INPUTS:
%   j       index of dynamics equation
%
% OUTPUTS:
%   Q       process noise covariance (6 x 6)
%
%==========================================================================
function Q = process_noise_covariance(j)
    Q = 0.1*eye(6);
end



%==========================================================================
% Measurement noise covariance.
%==========================================================================
%
% INPUTS:
%   i       index of measurement equation
%
% OUTPUTS:
%   R       measurement noise covariance (6 x 6)
%
%==========================================================================
function R = measurement_noise_covariance(i)
    R = eye(6);
end