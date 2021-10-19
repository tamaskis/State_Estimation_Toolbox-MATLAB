%% ps3
% Problem Set 3
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
pp = PLOT_PARAMETERS;

% random seed
seed = 2;



%% PROBLEM 2

% time parameters
dt = 1;         % time step [s]
tf = 10;        % simulation end time [s]
t = 0:dt:tf;    % time vector [s]
T = length(t);  % length of time vector

% system matrices
A = @(x,u,t) [eye(2),eye(2)*dt;zeros(2),eye(2)];
B = @(x,u,t) [zeros(2);eye(2)*dt];
C = @(x,t) [eye(2),zeros(2)];

% state (n) and measurement (m) dimensions
n = 4;
m = 2;

% process (Q) and measurement (R) noise covariances
Q = [zeros(2),zeros(2);zeros(2),eye(2)];
R = 9*eye(m);

% control input
u = [-2.5*cos(0.05*t);-2.5*sin(0.05*t)];
 
% initial prior distribution (i.e. initial state estimate and covariance)
mu0 = [1500;100;0;55];
Sigma0 = [250000*eye(2),zeros(2);zeros(2),eye(2)];

% ground truth initial condition
p0 = [1000;0];  % position [m]
s0 = [0;50];    % velocity [m/s]
x0 = [p0;s0];   % state vector

% initial condition structure for simulation
IC.x0_true = x0;

% ground truth simulation
[x,y] = simulate_linear(A,B,C,Q,R,u,IC,seed);

% state trajectory with measurement
figure('position',pp.plot_position);
hold on;
plot(x(1,:),x(2,:),'color',pp.cardinal_red,'linewidth',pp.line_width);
plot(y(1,2:end),y(2,2:end),'k:','linewidth',pp.line_width);
hold off;
grid on;
xlabel('$p_{1}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$p_{2}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
legend('true trajectory, $\mathbf{p}$',...
    'position measurement, $\mathbf{y}$','interpreter','latex',...
    'fontsize',pp.legend_font_size);



%% PROBLEM 3

% runs Kalman filter
[mu,Sigma] = KF(A,B,C,Q,R,t,u,y,mu0,Sigma0);

% plots estimated trajectory with position error ellipses (part (b))
figure('position',pp.plot_position);
hold on;
for i = 2:length(t)
    [x_ellipse,y_ellipse] = error_ellipse(mu(1:2,i),Sigma(1:2,1:2,i),0.95);
    if i == length(t)
        fill(x_ellipse,y_ellipse,'','facecolor',[0.75,0.75,0.75],...
            'edgecolor','none');
    else
        fill(x_ellipse,y_ellipse,'','facecolor',[0.75,0.75,0.75],...
            'edgecolor','none','handlevisibility','off');
    end
end
plot(x(1,:),x(2,:),'k:','linewidth',pp.line_width);
plot(mu(1,2:end),mu(2,2:end),'linewidth',pp.line_width,'color',...
    pp.cardinal_red);
hold off;
grid on;
xlabel('$p_{1}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$p_{2}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
legend('error ellipses for position estimate','true trajectory',...
    'position estimate','interpreter','latex','fontsize',...
    pp.legend_font_size,'location','best');


% plots true trajectory with velocity + velocity error ellipses (part (c))
figure('position',pp.plot_position);
hold on;
for i = 2:length(t)
    [x_ellipse,y_ellipse] = error_ellipse(mu(3:4,i),Sigma(3:4,3:4,i),0.95);
    x_ellipse = x_ellipse+x(1,i);
    y_ellipse = y_ellipse+x(2,i);
    if i == length(t)
        fill(x_ellipse,y_ellipse,'','facecolor',[0.75,0.75,0.75],...
            'edgecolor','none');
        plot(x(1,:),x(2,:),'k:','linewidth',pp.line_width);
        quiver(x(1,i),x(2,i),(10/9)*mu(3,i),(10/9)*mu(4,i),'linewidth',...
            pp.line_width,'color',pp.cardinal_red);
    else
        fill(x_ellipse,y_ellipse,'','facecolor',[0.75,0.75,0.75],...
            'edgecolor','none','handlevisibility','off');
        quiver(x(1,i),x(2,i),(10/9)*mu(3,i),(10/9)*mu(4,i),'linewidth',...
            pp.line_width,'color',pp.cardinal_red,'handlevisibility',...
            'off');
    end
end
hold off;
grid on;
xlabel('$p_{1}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$p_{2}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
legend('error ellipses for velocity estimate','true trajectory',...
    'velocity estimate','interpreter','latex','fontsize',...
    pp.legend_font_size,'location','best');



%% PROBLEM 4

% updated observation matrix
C = @(x,t) [zeros(2),eye(2)];

% updated initial prior distribution
mu0 = [1000;0;0;50];
Sigma0 = [eye(2),zeros(2);zeros(2),eye(2)];

% runs ground truth simulation + Kalman filter
[x,y] = simulate_linear(A,B,C,Q,R,u,IC,seed);
[mu,Sigma] = KF(A,B,C,Q,R,t,u,y,mu0,Sigma0);

% plots estimated trajectory with position error ellipses
figure('position',pp.plot_position);
hold on;
for i = 2:length(t)
    [x_ellipse,y_ellipse] = error_ellipse(mu(1:2,i),Sigma(1:2,1:2,i),0.95);
    if i == length(t)
        fill(x_ellipse,y_ellipse,'','facecolor',[0.75,0.75,0.75],...
            'edgecolor','none');
    else
        fill(x_ellipse,y_ellipse,'','facecolor',[0.75,0.75,0.75],...
            'edgecolor','none','handlevisibility','off');
    end
end
plot(x(1,:),x(2,:),'k:','linewidth',pp.line_width);
plot(mu(1,2:end),mu(2,2:end),'linewidth',pp.line_width,'color',...
    pp.cardinal_red);
hold off;
grid on;
xlabel('$p_{1}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$p_{2}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
legend('error ellipses for position estimate','true trajectory',...
    'position estimate','interpreter','latex','fontsize',...
    pp.legend_font_size,'location','best');

% plots true trajectory with velocity + velocity error ellipses
figure('position',pp.plot_position);
hold on;
for i = 2:length(t)
    [x_ellipse,y_ellipse] = error_ellipse(mu(3:4,i),Sigma(3:4,3:4,i),0.95);
    x_ellipse = x_ellipse+x(1,i);
    y_ellipse = y_ellipse+x(2,i);
    if i == length(t)
        fill(x_ellipse,y_ellipse,'','facecolor',[0.75,0.75,0.75],...
            'edgecolor','none');
        plot(x(1,:),x(2,:),'k:','linewidth',pp.line_width);
        quiver(x(1,i),x(2,i),(10/9)*mu(3,i),(10/9)*mu(4,i),'linewidth',...
            pp.line_width,'color',pp.cardinal_red);
    else
        fill(x_ellipse,y_ellipse,'','facecolor',[0.75,0.75,0.75],...
            'edgecolor','none','handlevisibility','off');
        quiver(x(1,i),x(2,i),(10/9)*mu(3,i),(10/9)*mu(4,i),'linewidth',...
            pp.line_width,'color',pp.cardinal_red,'handlevisibility',...
            'off');
    end
end
hold off;
grid on;
xlabel('$p_{1}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('$p_{2}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
legend('error ellipses for velocity estimate','true trajectory',...
    'velocity estimate','interpreter','latex','fontsize',...
    pp.legend_font_size,'location','best');