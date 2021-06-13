%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Problem Set 5 Problem 1

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

% roll mode parameters
a = 0.5;
b = 2;

% initial condition
%p0 = 1; % roll rate [rad/s]

% dimension parameters
n = 3;
k = 1;

% noise covariances
Q = [0.1,0,0;0,0,0;0,0,0];
R = 0.1;

% time parameters
dt = 0.1; % time step [s]
tf = 20; % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% simulated noise
w = [gaussian_random_sample(0,Q(1,1),T-1);zeros(1,T-1);zeros(1,T-1)];
v = gaussian_random_sample(zeros(k,1),R,T-1);

% control signal
u = sin(t);

% array preallocation
x = zeros(n,T);
y = zeros(k,T);

% samples initial condition for roll rate [rad/s]
p0 = gaussian_random_sample(0,0.1);

% stores initial condition
x(:,1) = [p0;a;b];

% functions for nonlinear dynamics and nonlinear measurements
f = @(x,u) f_discrete(x,u,dt);
g = @(x) g_discrete(x);

% simulation
for tt = 1:(T-1)
    
    % propagates state vector with noise
    x(:,tt+1) = f(x(:,tt),u(:,tt))+w(:,tt);
    
    % measurement with noise
    y(:,tt+1) = g(x(:,tt))+v(:,tt);
    
end

% initializes figure for roll rate
figure('position',two_subplot_position);

% true roll rate
subplot(1,2,1);
plot(t,x(1,:),'color',cardinal_red,'linewidth',line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$p\;[\mathrm{rad/s}]$','interpreter','latex','fontsize',...
    axis_font_size);
title('\textbf{True Roll Rate}','interpreter','latex','fontsize',...
    title_font_size)

% measured roll rate
subplot(1,2,2);
plot(t(2:end),y(1,2:end),'color',cardinal_red,'linewidth',line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$p\;[\mathrm{rad/s}]$','interpreter','latex','fontsize',...
    axis_font_size);
title('\textbf{Measured Roll Rate}','interpreter','latex','fontsize',...
    title_font_size)



%% FILTERING

% initial state estimate
mu0 = [0;0;0];
Sigma0 = 10*eye(3);

% dynamics and measurement Jacobians
A = @(x,u) dynamics_jacobian(x,u,dt);
C = @(x) measurement_jacobian(x);

% redefines Q for EKF
Q = [0.1,0,0;0,0.01,0;0,0,0.01];

% runs EKF
[mu,Sigma] = EKF(f,g,A,C,Q,R,u,y,mu0,Sigma0);

% 95% confidence interval
mu_plus = zeros(size(mu));
mu_minus = zeros(size(mu));
for i = 1:tt
    mu_plus(:,i) = mu(:,i)+1.96*sqrt(diag(Sigma(:,:,i)));
    mu_minus(:,i) = mu(:,i)-1.96*sqrt(diag(Sigma(:,:,i)));
end

% initializes figure
figure('position',three_subplot_position);

% true roll rate and its estimate
subplot(1,3,1);
hold on;
patch([t,fliplr(t)],[mu_plus(1,:),fliplr(mu_minus(1,:))],...
    matlab_light_red,'edgecolor','none','handlevisibility','off');
plot(t,x(1,:),'linewidth',line_width);
plot(t(2:end),mu(1,2:end),'linewidth',line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$p\;[\mathrm{rad/s}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true roll rate','estimated roll rate (with 95\% C.I.)',...
    'interpreter','latex','fontsize',legend_font_size);

% true a and its estimate
subplot(1,3,2);
hold on;
patch([t,fliplr(t)],[mu_plus(2,:),fliplr(mu_minus(2,:))],...
    matlab_light_red,'edgecolor','none','handlevisibility','off');
plot(t,x(2,:),'linewidth',line_width);
plot(t(2:end),mu(2,2:end),'linewidth',line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$a$','interpreter','latex','fontsize',axis_font_size);
legend('true $a$','estimated $a$ (with 95\% C.I.)','interpreter',...
    'latex','fontsize',legend_font_size);

% true b and its estimate
subplot(1,3,3);
hold on;
patch([t,fliplr(t)],[mu_plus(3,:),fliplr(mu_minus(3,:))],...
    matlab_light_red,'edgecolor','none','handlevisibility','off');
plot(t,x(3,:),'linewidth',line_width);
plot(t(2:end),mu(3,2:end),'linewidth',line_width);
hold off;
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$b$','interpreter','latex','fontsize',axis_font_size);
legend('true $b$','estimated $b$ (with 95\% C.I.)','interpreter',...
    'latex','fontsize',legend_font_size);



%% ADDITIONAL FUNCTIONS

%=========================================================================%
% Discrete-time nonlinear dynamics.
%=========================================================================%
% INPUT: x - state vector at current time step
%        u - control input at current time step
%        dt - time step [s]
% OUTPUT: f_eval - evaluation of f (produces state vector at next time 
%                  step)
function f_eval = f_discrete(x,u,dt)

    % unpacks state vector
    p = x(1);
    a = x(2);
    b = x(3);
    
    % evalutes f(xt,ut)
    f_eval = [p+(-a*p+b*u)*dt;a;b];
    
end



%=========================================================================%
% Discrete-time nonlinear measurement.
%=========================================================================%
% INPUT: x - state vector at current time step
% OUTPUT: g_eval - evaluation of g (produces measurement at current time 
%                  step)
function g_eval = g_discrete(x)

    % extracts "p" from state vector
    p = x(1);
    
    % evalutes g(xt,ut)
    g_eval = p;
    
end



%=========================================================================%
% Dynamics Jacobian.
%=========================================================================%
% INPUT: x - state vector
%        u - control input
%        dt - time step [s]
% OUTPUT: A - dynamics Jacobian
function A = dynamics_jacobian(x,u,dt)

    % extacts "p" and "a" from state vector
    p = x(1);
    a = x(2);
    
    % assembles dynamics Jacobian
    A = [1-a*dt   -p*dt   u*dt;
         0         1      0;
         0         0      1];
    
end



%=========================================================================%
% Measurement Jacobian.
%=========================================================================%
% INPUT: x - state vector
% OUTPUT: C - measurement Jacobian
function C = measurement_jacobian(x)
    C = [1,0,0];
end