%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Problem Set 5 Problem 3

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

% time parameters
dt = 0.1; % time step [s]
tf = 20; % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% dimension parameters
n = 3;
k = 1;

% noise covariances
Q = 0.1*eye(3)*dt;
R = 0.1;

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
f = @(x,u) f_discrete(x,u,dt);
g = @(x) g_discrete(x);

% simulation
for tt = 1:(T-1)
    
    % propagates state vector with noise
    x(:,tt+1) = f(x(:,tt),u(:,tt))+w(:,tt);
    
    % measurement with noise
    y(:,tt+1) = g(x(:,tt))+v(:,tt);
    
end

% initializes figure for dynamics simulation
figure('position',three_subplot_position);

% x-position
subplot(1,3,1);
plot(t,x(1,:),'color',cardinal_red,'linewidth',line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$p_{x}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);

% y-position
subplot(1,3,2);
plot(t,x(2,:),'color',cardinal_red,'linewidth',line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$p_{y}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);

% heading angle
subplot(1,3,3);
plot(t,x(3,:),'color',cardinal_red,'linewidth',line_width);
grid on;
axis square;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\theta\;[\mathrm{rad}]$','interpreter','latex','fontsize',...
    axis_font_size);

% measurement
figure('position',plot_position)
plot(t(2:end),y(1,2:end),'color',cardinal_red,'linewidth',line_width);
grid on;
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('${y}_{t}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);



%% PART (A) - EKF

% dynamics and measurement Jacobians
A = @(x,u) dynamics_jacobian(x,u,dt);
C = @(x) measurement_jacobian(x);

% runs EKF
tic;
[mu,Sigma,Ob] = EKF(f,g,A,C,Q,R,u,y,mu0,Sigma0);
t_EKF = toc;

% finds rank of observability matrix at each time step
Ob_rank = zeros(size(t));
for i = 2:length(t)
    Ob_rank(i) = rank(Ob(:,:,i));
end

% 95% confidence interval
mu_plus = zeros(size(x));
mu_minus = zeros(size(x));
for i = 1:tt
    mu_plus(:,i) = mu(:,i)+1.96*sqrt(diag(Sigma(:,:,i)));
    mu_minus(:,i) = mu(:,i)-1.96*sqrt(diag(Sigma(:,:,i)));
end

% initializes figure
figure('position',three_subplot_position);

% true x-position and its estimate
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
ylabel('$p_{x}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true $x$-position','estimated $x$-position (with 95\% C.I.)',...
    'interpreter','latex','fontsize',legend_font_size);

% true y-position and its estimate
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
ylabel('$p_{y}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true $y$-position','estimated $y$-position (with 95\% C.I.)',...
    'interpreter','latex','fontsize',legend_font_size);

% true heading angle and its estimate
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
ylabel('$\theta\;[\mathrm{rad}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true heading angle','estimated heading angle (with 95\% C.I.)',...
    'interpreter',...
    'latex','fontsize',legend_font_size);

% rank of observability matrix
figure('position',plot_position);
plot(t(2:end),Ob_rank(2:end),'color',cardinal_red,'linewidth',line_width);
grid on;
ylim([0,3]);
xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$\mathrm{rank}(\mathcal{O})$','interpreter','latex','fontsize',...
    axis_font_size);



%% PART (B) - iEKF

% runs iEKF
tic;
[mu,Sigma] = iEKF(f,g,A,C,Q,R,u,y,mu0,Sigma0);
t_iEKF = toc;

% 95% confidence interval
mu_plus = zeros(size(x));
mu_minus = zeros(size(x));
for i = 1:tt
    mu_plus(:,i) = mu(:,i)+1.96*sqrt(diag(Sigma(:,:,i)));
    mu_minus(:,i) = mu(:,i)-1.96*sqrt(diag(Sigma(:,:,i)));
end

% initializes figure
figure('position',three_subplot_position);

% true x-position and its estimate
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
ylabel('$p_{x}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true $x$-position','estimated $x$-position (with 95\% C.I.)',...
    'interpreter','latex','fontsize',legend_font_size);

% true y-position and its estimate
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
ylabel('$p_{y}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true $y$-position','estimated $y$-position (with 95\% C.I.)',...
    'interpreter','latex','fontsize',legend_font_size);

% true heading angle and its estimate
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
ylabel('$\theta\;[\mathrm{rad}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true heading angle','estimated heading angle (with 95\% C.I.)',...
    'interpreter',...
    'latex','fontsize',legend_font_size);



%% PART (C) - UKF

% runs UKF
tic;
[mu,Sigma] = UKF(f,g,Q,R,u,y,mu0,Sigma0);
t_UKF = toc;

% 95% confidence interval
mu_plus = zeros(size(x));
mu_minus = zeros(size(x));
for i = 1:tt
    mu_plus(:,i) = mu(:,i)+1.96*sqrt(diag(Sigma(:,:,i)));
    mu_minus(:,i) = mu(:,i)-1.96*sqrt(diag(Sigma(:,:,i)));
end

% initializes figure
figure('position',three_subplot_position);

% true x-position and its estimate
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
ylabel('$p_{x}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true $x$-position','estimated $x$-position (with 95\% C.I.)',...
    'interpreter','latex','fontsize',legend_font_size);

% true y-position and its estimate
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
ylabel('$p_{y}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true $y$-position','estimated $y$-position (with 95\% C.I.)',...
    'interpreter','latex','fontsize',legend_font_size);

% true heading angle and its estimate
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
ylabel('$\theta\;[\mathrm{rad}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true heading angle','estimated heading angle (with 95\% C.I.)',...
    'interpreter',...
    'latex','fontsize',legend_font_size);



%% PART (D) SEQUENTIAL IMPORTANCE RESAMPLING PARTICLE FILTER

% number of particles
N = 1000;

% produces initial particles (random samples of prior distribution)
x0 = gaussian_random_sample(mu0,Sigma0,N);

% function handle for sampling process noise
w = @() gaussian_random_sample(zeros(n,1),Q);

% function handle for measurement noise PDF
pv = @(x) mvnpdf(x,zeros(k,1),R);

% runs SIR particle filter
tic;
[mu,Sigma] = SIR(f,g,w,pv,u,y,x0,N);
t_SIR = toc;

% 95% confidence interval
mu_plus = zeros(size(mu));
mu_minus = zeros(size(mu));
for i = 1:tt
    mu_plus(:,i) = mu(:,i)+1.96*sqrt(diag(Sigma(:,:,i)));
    mu_minus(:,i) = mu(:,i)-1.96*sqrt(diag(Sigma(:,:,i)));
end

% initializes figure
figure('position',three_subplot_position);

% true x-position and its estimate
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
ylabel('$p_{x}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true $x$-position','estimated $x$-position (with 95\% C.I.)',...
    'interpreter','latex','fontsize',legend_font_size);

% true y-position and its estimate
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
ylabel('$p_{y}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true $y$-position','estimated $y$-position (with 95\% C.I.)',...
    'interpreter','latex','fontsize',legend_font_size);

% true heading angle and its estimate
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
ylabel('$\theta\;[\mathrm{rad}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true heading angle','estimated heading angle (with 95\% C.I.)',...
    'interpreter',...
    'latex','fontsize',legend_font_size);



%% PART (F) AVERAGE COMPUTATION TIME

% converts measured times to averages based on T (number of iterations)
t_EKF = t_EKF/T;
t_iEKF = t_iEKF/T;
t_UKF = t_UKF/T;
t_SIR = t_SIR/T;

% prints results
format shortE;
table(t_EKF,t_iEKF,t_UKF,t_SIR)
format short;




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
% INPUT: x - state vector at current time step
% OUTPUT: g_eval - evaluation of g (produces measurement at current time 
%                  step)
function g_eval = g_discrete(x)

    % extracts "p" from state vector
    p = x(1:2);
    
    % evalutes g(xt,ut)
    g_eval = norm(p);
    
end



%=========================================================================%
% Dynamics Jacobian.
%=========================================================================%
% INPUT: x - state vector
%        u - control input
%        dt - time step [s]
% OUTPUT: A - dynamics Jacobian
function A = dynamics_jacobian(x,u,dt)

    % extacts "theta" from state vector
    theta = x(3);
    
    % extacts "nu" from control input
    nu = u(1);
    
    % assembles dynamics Jacobian
    A = [1   0   -nu*sin(theta)*dt;
         0   1    nu*cos(theta)*dt;
         0   0    1];
    
end



%=========================================================================%
% Measurement Jacobian.
%=========================================================================%
% INPUT: x - state vector
% OUTPUT: C - measurement Jacobian
function C = measurement_jacobian(x)
    
    % extacts "px" and "py" from state vector
    px = x(1);
    py = x(2);
    
    % assembles measurement Jacobian
    C = [px/sqrt(px^2+py^2)   py/sqrt(px^2+py^2)   0];
    
end