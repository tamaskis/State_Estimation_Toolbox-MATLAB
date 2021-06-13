%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Problem Set 6 Problems 3 and 4

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
tf = 100; % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% dimension parameters
n = 11;
k = 8;

% noise covariances
Q = [0.1*eye(3)*dt,zeros(3,8);zeros(8,3),0.01*eye(8)];
R = 0.1*eye(k);

% true feature locations
m1_true = [0;0];
m2_true = [10;0];
m3_true = [10;10];
m4_true = [0;10];

% initial state estimate (prior distribution)
mu0 = [0;0;0;m1_true+1;m2_true+1;m3_true+1;m4_true+1];
Sigma0 = 0.001*eye(length(mu0));

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
g = @(x,t) g_discrete(x,t);

% simulation (using true feature locations)
for tt = 1:(T-1)
    
    % propagates state vector with noise
    x(:,tt+1) = f([x(1:3,tt);m1_true;m2_true;m3_true;m4_true],u(:,tt),...
        tt)+w(:,tt);
    
    % measurement with noise
    y(:,tt) = g([x(1:3,tt);m1_true;m2_true;m3_true;m4_true],tt)+v(:,tt);
    
end

% extracts time history of true position
px_true = x(1,:);
py_true = x(2,:);

% dynamics and measurement Jacobians
A = @(x,u,t) dynamics_jacobian(x,u,t,dt);
C = @(x,t) measurement_jacobian(x,t);



%% PROBLEM 3(A)

% EKF SLAM
[mu,Sigma] = EKF_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,5:8);

% extracts time history of position estimates
px = mu(1,2:end);
py = mu(2,2:end);

% extracts time history of feature location estimates
m1 = mu(4:5,2:end);
m2 = mu(6:7,2:end);
m3 = mu(8:9,2:end);
m4 = mu(10:11,2:end);

% map
figure('position',plot_position);
hold on;
for i = 1:5:length(t)
    [x_ell,y_ell] = error_ellipse(mu(1:2,i),Sigma(1:2,1:2,i),0.95);
    patch(x_ell,y_ell,matlab_light_red,'edgecolor','none',...
        'handlevisibility','off');
end
plot(px_true,py_true,'k:','linewidth',line_width);
plot(px,py,'linewidth',line_width,'color',cardinal_red);
scatter([m1_true(1),m2_true(1),m3_true(1),m4_true(1)],[m1_true(2),...
    m2_true(2),m3_true(2),m4_true(2)],50,...
    'markeredgecolor','k','MarkerFaceColor','g','linewidth',1);
scatter([mu0(4),mu0(6),mu0(8),mu0(10)],[mu0(5),mu0(7),mu0(9),mu0(11)],...
    50,'markeredgecolor','k','MarkerFaceColor','r','linewidth',1);
scatter(m1(1,:),m1(2,:),5,matlab_blue,'filled');
scatter(m2(1,:),m2(2,:),5,matlab_red,'filled');
scatter(m3(1,:),m3(2,:),5,matlab_yellow,'filled');
scatter(m4(1,:),m4(2,:),5,matlab_purple,'filled');
hold off;
grid on;
axis equal;
xlim([-5,30]);
ylim([-5,30]);
xlabel('$p_{x}$','interpreter','latex','fontsize',axis_font_size);
ylabel('$p_{y}$','interpreter','latex','fontsize',axis_font_size);
legend('true trajectory','estimated trajectory (with 95\% C.I.)',...
    'true feature locations','initial guess of feature locations',...
    'feature 1 estimate','feature 2 estimate','feature 3 estimate',...
    'feature 4 estimate','interpreter','latex','fontsize',...
    legend_font_size,'location','northeast');



%% PROBLEM 3(B)

% iEKF SLAM
[mu,Sigma] = iEKF_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,5:8);

% extracts time history of position estimates
px = mu(1,2:end);
py = mu(2,2:end);

% extracts time history of feature location estimates
m1 = mu(4:5,2:end);
m2 = mu(6:7,2:end);
m3 = mu(8:9,2:end);
m4 = mu(10:11,2:end);

% map
figure('position',plot_position);
hold on;
for i = 1:5:length(t)
    [x_ell,y_ell] = error_ellipse(mu(1:2,i),Sigma(1:2,1:2,i),0.95);
    patch(x_ell,y_ell,matlab_light_red,'edgecolor','none',...
        'handlevisibility','off');
end
plot(px_true,py_true,'k:','linewidth',line_width);
plot(px,py,'linewidth',line_width,'color',cardinal_red);
scatter([m1_true(1),m2_true(1),m3_true(1),m4_true(1)],[m1_true(2),...
    m2_true(2),m3_true(2),m4_true(2)],50,...
    'markeredgecolor','k','MarkerFaceColor','g','linewidth',1);
scatter([mu0(4),mu0(6),mu0(8),mu0(10)],[mu0(5),mu0(7),mu0(9),mu0(11)],...
    50,'markeredgecolor','k','MarkerFaceColor','r','linewidth',1);
scatter(m1(1,:),m1(2,:),5,matlab_blue,'filled');
scatter(m2(1,:),m2(2,:),5,matlab_red,'filled');
scatter(m3(1,:),m3(2,:),5,matlab_yellow,'filled');
scatter(m4(1,:),m4(2,:),5,matlab_purple,'filled');
hold off;
grid on;
axis equal;
xlim([-5,30]);
ylim([-5,30]);
xlabel('$p_{x}$','interpreter','latex','fontsize',axis_font_size);
ylabel('$p_{y}$','interpreter','latex','fontsize',axis_font_size);
legend('true trajectory','estimated trajectory (with 95\% C.I.)',...
    'true feature locations','initial guess of feature locations',...
    'feature 1 estimate','feature 2 estimate','feature 3 estimate',...
    'feature 4 estimate','interpreter','latex','fontsize',...
    legend_font_size,'location','northeast');



%% PROBLEM 4(A)

% EKS SLAM
[mu,Sigma] = EKS_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,5:8);

% extracts time history of position estimates
px = mu(1,2:end);
py = mu(2,2:end);

% extracts time history of feature location estimates
m1 = mu(4:5,2:end);
m2 = mu(6:7,2:end);
m3 = mu(8:9,2:end);
m4 = mu(10:11,2:end);

% +/- 2-sigma bounds
[lower_bound,upper_bound] = sigma_bounds(mu,Sigma,2);

% upper and lower 2-sigma bounds for position estimates
px_minus = lower_bound(1,:);
px_plus = upper_bound(1,:);
py_minus = lower_bound(2,:);
py_plus = upper_bound(2,:);

% map
figure('position',plot_position);
hold on;
for i = 1:5:length(t)
    [x_ell,y_ell] = error_ellipse(mu(1:2,i),Sigma(1:2,1:2,i),0.95);
    patch(x_ell,y_ell,matlab_light_red,'edgecolor','none',...
        'handlevisibility','off');
end
plot(px_true,py_true,'k:','linewidth',line_width);
plot(px,py,'linewidth',line_width,'color',cardinal_red);
scatter([m1_true(1),m2_true(1),m3_true(1),m4_true(1)],[m1_true(2),...
    m2_true(2),m3_true(2),m4_true(2)],50,...
    'markeredgecolor','k','MarkerFaceColor','g','linewidth',1);
scatter([mu0(4),mu0(6),mu0(8),mu0(10)],[mu0(5),mu0(7),mu0(9),mu0(11)],...
    50,'markeredgecolor','k','MarkerFaceColor','r','linewidth',1);
scatter(m1(1,:),m1(2,:),5,matlab_blue,'filled');
scatter(m2(1,:),m2(2,:),5,matlab_red,'filled');
scatter(m3(1,:),m3(2,:),5,matlab_yellow,'filled');
scatter(m4(1,:),m4(2,:),5,matlab_purple,'filled');
hold off;
grid on;
axis equal;
xlim([-5,30]);
ylim([-5,30]);
xlabel('$p_{x}$','interpreter','latex','fontsize',axis_font_size);
ylabel('$p_{y}$','interpreter','latex','fontsize',axis_font_size);
legend('true trajectory','estimated trajectory (with 95\% C.I.)',...
    'true feature locations','initial guess of feature locations',...
    'feature 1 estimate','feature 2 estimate','feature 3 estimate',...
    'feature 4 estimate','interpreter','latex','fontsize',...
    legend_font_size,'location','northeast');



%% PROBLEM 4(B)

% iEKS SLAM
[mu,Sigma] = iEKS_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,5:8);

% extracts time history of position estimates
px = mu(1,2:end);
py = mu(2,2:end);

% extracts time history of feature location estimates
m1 = mu(4:5,2:end);
m2 = mu(6:7,2:end);
m3 = mu(8:9,2:end);
m4 = mu(10:11,2:end);

% map
figure('position',plot_position);
hold on;
for i = 1:5:length(t)
    [x_ell,y_ell] = error_ellipse(mu(1:2,i),Sigma(1:2,1:2,i),0.95);
    patch(x_ell,y_ell,matlab_light_red,'edgecolor','none',...
        'handlevisibility','off');
end
plot(px_true,py_true,'k:','linewidth',line_width);
plot(px,py,'linewidth',line_width,'color',cardinal_red);
scatter([m1_true(1),m2_true(1),m3_true(1),m4_true(1)],[m1_true(2),...
    m2_true(2),m3_true(2),m4_true(2)],50,...
    'markeredgecolor','k','MarkerFaceColor','g','linewidth',1);
scatter([mu0(4),mu0(6),mu0(8),mu0(10)],[mu0(5),mu0(7),mu0(9),mu0(11)],...
    50,'markeredgecolor','k','MarkerFaceColor','r','linewidth',1);
scatter(m1(1,:),m1(2,:),5,matlab_blue,'filled');
scatter(m2(1,:),m2(2,:),5,matlab_red,'filled');
scatter(m3(1,:),m3(2,:),5,matlab_yellow,'filled');
scatter(m4(1,:),m4(2,:),5,matlab_purple,'filled');
hold off;
grid on;
axis equal;
xlim([-5,30]);
ylim([-5,30]);
xlabel('$p_{x}$','interpreter','latex','fontsize',axis_font_size);
ylabel('$p_{y}$','interpreter','latex','fontsize',axis_font_size);
legend('true trajectory','estimated trajectory (with 95\% C.I.)',...
    'true feature locations','initial guess of feature locations',...
    'feature 1 estimate','feature 2 estimate','feature 3 estimate',...
    'feature 4 estimate','interpreter','latex','fontsize',...
    legend_font_size,'location','northeast');



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
    m1x = x(4);
    m1y = x(5);
    m2x = x(6);
    m2y = x(7);
    m3x = x(8);
    m3y = x(9);
    m4x = x(10);
    m4y = x(11);
    
    % unpacks control vector
    nu = u(1);
    phi = u(2);
    
    % evalutes f(xt,ut)
    f_eval = [px+nu*cos(theta)*dt;
              py+nu*sin(theta)*dt;
              theta+phi*dt;
              m1x;
              m1y;
              m2x;
              m2y;
              m3x;
              m3y;
              m4x;
              m4y];
    
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
function g_eval = g_discrete(x,t)
    
    % unpacks state vector
    px = x(1);
    py = x(2);
    theta = x(3);
    m1x = x(4);
    m1y = x(5);
    m2x = x(6);
    m2y = x(7);
    m3x = x(8);
    m3y = x(9);
    m4x = x(10);
    m4y = x(11);
    
    % evalutes g(xt,ut)
    g_eval = [sqrt((m1x-px)^2+(m1y-py)^2);
              sqrt((m2x-px)^2+(m2y-py)^2);
              sqrt((m3x-px)^2+(m3y-py)^2);
              sqrt((m4x-px)^2+(m4y-py)^2);
              atan2(m1y-py,m1x-px)-theta;
              atan2(m2y-py,m2x-px)-theta;
              atan2(m3y-py,m3x-px)-theta;
              atan2(m4y-py,m4x-px)-theta];
                      
end



%=========================================================================%
% Dynamics Jacobian.
%=========================================================================%
%
% INPUTS:
%   x       state vector (11 x 1)
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
    
    % problem 1a dynamics Jacobian
    A1a = [1   0   -nu*sin(theta)*dt;
           0   1    nu*cos(theta)*dt;
           0   0    1];
    
    % problem 2 dynamics Jacobian
    A2 = eye(8);
    
    % assembles full dynamics Jacobian
    A = [A1a          zeros(3,8);
         zeros(8,3)   A2];
    
end



%=========================================================================%
% Measurement Jacobian.
%=========================================================================%
%
% INPUTS:
%   x       state vector (11 x 1)
%   t       current iteration (corresponding to discrete time)
%
% OUTPUTS:
%   C       measurement Jacobian
%
%=========================================================================%
function C = measurement_jacobian(x,t)
    
    % creates "mx" and "my" vectors from state vector
    mx = x([4,6,8,9]);
    my = x([5,7,9,11]);
    
    % unpacks state vector
    px = x(1);
    py = x(2);
    m1x = x(4);
    m1y = x(5);
    m2x = x(6);
    m2y = x(7);
    m3x = x(8);
    m3y = x(9);
    m4x = x(10);
    m4y = x(11);
    
    % problem 1a measurement Jacobian
    C1a = [(px-m1x)/sqrt((m1x-px)^2+(m1y-py)^2)   (py-m1y)/sqrt((m1x-px)^2+(m1y-py)^2)   0;
           (px-m2x)/sqrt((m2x-px)^2+(m2y-py)^2)   (py-m2y)/sqrt((m2x-px)^2+(m2y-py)^2)   0;
           (px-m3x)/sqrt((m3x-px)^2+(m3y-py)^2)   (py-m3y)/sqrt((m3x-px)^2+(m3y-py)^2)   0;
           (px-m4x)/sqrt((m4x-px)^2+(m4y-py)^2)   (py-m4y)/sqrt((m4x-px)^2+(m4y-py)^2)   0];
    
    % problem 1b measurement Jacobian
    C1b = [(m1y-py)/((m1x-px)^2+(m1y-py)^2)   (px-m1x)/((m1x-px)^2+(m1y-py)^2)   -1;
           (m2y-py)/((m2x-px)^2+(m2y-py)^2)   (px-m2x)/((m2x-px)^2+(m2y-py)^2)   -1;
           (m3y-py)/((m3x-px)^2+(m3y-py)^2)   (px-m3x)/((m3x-px)^2+(m3y-py)^2)   -1;
           (m4y-py)/((m4x-px)^2+(m4y-py)^2)   (px-m4x)/((m4x-px)^2+(m4y-py)^2)   -1];
    
    % problem 2 measurement Jacobian
    C2 = zeros(4,8);
    for i = 1:4
        for j = 1:4
            if i == j
                C2(i,2*j-1) = (py-my(j))/((mx(j)-px)^2+(my(j)-py)^2);
                C2(i,2*j) = (mx(j)-px)/((mx(j)-px)^2+(my(j)-py)^2);
            end
        end
    end
    
    % problem 3 measurement Jacobian
    C3 = zeros(4,8);
    for i = 1:4
        for j = 1:4
            if i == j
                C3(i,2*j-1) = (mx(j)-px)/sqrt((mx(j)-px)^2+(my(j)-py)^2);
                C3(i,2*j) = (my(j)-py)/sqrt((mx(j)-px)^2+(my(j)-py)^2);
            end
        end
    end
    
    % assembles full measurement Jacobian
    C = [C1a   C3;
         C1b   C2];
    
end