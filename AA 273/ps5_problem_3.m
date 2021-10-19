%% ps5_problem_3
% Problem Set 5, Problem 3
% AA 273 - State Estimation and Filtering for Robotic Perception
%
% Author: Tamas Kis
% Last Update: 2021-08-18



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear;clc;close all;

% adds path to root directory and all subdirectories
addpath(genpath("../"));

% loads plot parameters
pp = PLOT_PARAMETERS();

% random seed
seed = 2;



%% SIMULATION

% time parameters
dt = 0.1;   % time step [s]
tf = 20;    % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% dimension parameters
n = 3;
k = 1;

% noise covariances
Q = 0.1*eye(n)*dt;
R = 0.1;

% initial state estimate (prior distribution)
mu0 = [0;0;0];
Sigma0 = 0.01*eye(n);

% control signal
u = [ones(size(t));sin(t)];

% functions for nonlinear dynamics and nonlinear measurements
f = @(x,u,t) discrete_dynamics(x,u,t,dt);
g = @(x,t) discrete_measurement(x,t);

% initial condition structure for simulation
IC.x0 = mu0;
IC.P0 = Sigma0;

% ground truth simulation
[x,y] = simulate_nonlinear(f,g,Q,R,t,u,IC,seed);

% plots simulation results
plot_simulation_results(pp,t,x,y)



%% PART (A) - EKF

% dynamics and measurement Jacobians
A = @(x,u,t) dynamics_jacobian(x,u,t,dt);
C = @(x,t) measurement_jacobian(x,t);

% runs EKF
[mu,Sigma,t_EKF,rank_Ob] = EKF(f,g,A,C,Q,R,t,u,y,mu0,Sigma0);

% plots results
plot_filter_results(pp,"EKF",t,x,mu,Sigma,rank_Ob);



%% PART (B) - iEKF

% runs iEKF
[mu,Sigma,t_iEKF,rank_Ob] = EKF(f,g,A,C,Q,R,t,u,y,mu0,Sigma0);

% plots results
plot_filter_results(pp,"iEKF",t,x,mu,Sigma,rank_Ob);



%% PART (C) - UKF

% runs UKF
[mu,Sigma,t_UKF] = UKF(f,g,Q,R,t,u,y,mu0,Sigma0);

% plots results
plot_filter_results(pp,"UKF",t,x,mu,Sigma);



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
[mu,Sigma] = SIR(f,g,w,pv,t,u,y,x0,N);
t_SIR = toc;

% plots results
plot_filter_results(pp,"Particle Filter (SIR)",t,x,mu,Sigma);



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

%==========================================================================
% discrete_dynamics  Discrete dynamics equation.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (3×1 double) state vector at current sample time
%   u       - (2×1 double) control input at current sample time
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
    px = x(1);
    py = x(2);
    theta = x(3);
    
    % unpacks control vector
    nu = u(1);
    phi = u(2);
    
    % evalutes f(x,u,t)
    x_next = [px+nu*cos(theta)*dt;py+nu*sin(theta)*dt;theta+phi*dt];
    
end



%==========================================================================
% discrete_measurement  Discrete measurement equation.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x 	- (3×1 double) state vector
%   t  	- (1×1 double) current time [s]
%
% -------
% OUTPUT:
% -------
%   y 	- (1×1 double) measurement
%
%==========================================================================
function y = discrete_measurement(x,t)

    % extracts "p" from state vector
    p = x(1:2);
    
    % evalutes g(xt,ut)
    y = norm(p);
    
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
%   t   - (1×1 double) time [s]
%   dt  - (1×1 double) time step [s]
%
% -------
% OUTPUT:
% -------
%   A   - (3×3 double) dynamics Jacobian
%
%==========================================================================
function A = dynamics_jacobian(x,u,t,dt)

    % extacts "theta" from state vector
    theta = x(3);
    
    % extacts "nu" from control input
    nu = u(1);
    
    % assembles dynamics Jacobian
    A = [1   0   -nu*sin(theta)*dt;
         0   1    nu*cos(theta)*dt;
         0   0    1];
    
end



%==========================================================================
% measurement_jacobian  Measurement Jacobian.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x   - (3×1 double) state vector
%   t  	- (1×1 double) time [s]
%
% -------
% OUTPUT:
% -------
%   C   - (1×3) measurement Jacobian
%
%==========================================================================
function C = measurement_jacobian(x,t)
    
    % extacts "px" and "py" from state vector
    px = x(1);
    py = x(2);
    
    % assembles measurement Jacobian
    C = [px/sqrt(px^2+py^2)   py/sqrt(px^2+py^2)   0];
    
end



%==========================================================================
% plot_simulation_results  Plots the simulation results.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   pp      - (struct) structure storing plot parameters
%   t     	- (1×T double) time vector [s]
%   x       - (n×T double) ground truth state time history
%   y       - (m×T double) ground truth measurement time history
%
%==========================================================================
function plot_simulation_results(plot_parameters,t,x,y)
    
    % extracts plot parameters
    three_subplot_position = plot_parameters.three_subplot_position;
    plot_position = plot_parameters.plot_position;
    line_width = plot_parameters.line_width;
    axis_font_size = plot_parameters.axis_font_size;
    cardinal_red = plot_parameters.cardinal_red;
    
    % ----------------------------
    % Plot of dynamics simulation.
    % ----------------------------
    
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

    % -------------------------------
    % Plot of measurement simulation.
    % -------------------------------
    
    % measurement
    figure('position',plot_position)
    plot(t(2:end),y(1,2:end),'color',cardinal_red,'linewidth',line_width);
    grid on;
    xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
        axis_font_size);
    ylabel('${y}_{t}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
        axis_font_size);
    
end



%==========================================================================
% plot_filter_results  Plots the filter results.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   plot_parameters - (struct) structure storing plot parameters
%   filter          - (string) filter name
%   t           	- (1×T) time vector [s]
%   x           	- (n×T) ground truth state time history
%   mu           	- (n×T) state estimate time history
%   Sigma          	- (n×n×T) error covariance time history
%   rank_Ob     	- (T×1) (OPTIONAL) rank of the observability matrix
%
%==========================================================================
function plot_filter_results(plot_parameters,filter,t,x,mu,Sigma,rank_Ob)
    
    % extracts plot parameters
    three_subplot_position = plot_parameters.three_subplot_position;
    plot_position = plot_parameters.plot_position;
    line_width = plot_parameters.line_width;
    title_font_size = plot_parameters.title_font_size;
    axis_font_size = plot_parameters.axis_font_size;
    legend_font_size = plot_parameters.legend_font_size;
    matlab_light_red = plot_parameters.matlab_light_red;
    
    % 2-sigma bounds
    [mu_minus,mu_plus] = sigma_bounds(mu,Sigma,2);

    % ---------------------------------------
    % Plot of state estimate with true state.
    % ---------------------------------------
    
    % initializes figure
    figure('position',three_subplot_position);
    sgtitle("\textbf{"+filter+" Results}",'interpreter','latex',...
        'fontsize',title_font_size);

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
    legend('true $x$-position',...
        'estimated $x$-position (with $\pm2\sigma$ bounds)',...
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
    legend('true $y$-position',...
        'estimated $y$-position (with $\pm2\sigma$ bounds)',...
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
    legend('true heading angle',...
        'estimated heading angle (with $\pm2\sigma$ bounds)',...
        'interpreter','latex','fontsize',legend_font_size);

    % rank of observability matrix
    if nargin == 10
        figure('position',plot_position);
        plot(t(2:end),rank_Ob(2:end),'color',cardinal_red,'linewidth',...
            line_width);
        grid on;
        ylim([0,3]);
        xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
            axis_font_size);
        ylabel('$\mathrm{rank}(\mathcal{O})$','interpreter','latex',...
            'fontsize',axis_font_size);
        title("\textbf{"+filter+" Observability Analysis}",...
            'interpreter','latex','fontsize',title_font_size);
    end
    
end