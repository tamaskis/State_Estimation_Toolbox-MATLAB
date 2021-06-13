%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Problem Set 3

%-------------------------------------------------------------------------%



%% SCRIPT SETUP

% clears variables and command window, closes all figures
clear;
clc;
close all;

% adds path to Estimation Toolbox
addpath("..");

% plot parameters
plot_position = [540,300,700,500]; % plot position [x,y,l,w]
line_width = 1.5; % line width [#]
axis_font_size = 18; % axis label font size [#]
legend_font_size = 14; % legend font size [#]
cardinal_red = [140,21,21]/255; % color for plots [rgb]



%% PROBLEM 2(b)

% initial condition
p0 = [1000;0]; % position [m]
s0 = [0;50]; % velocity [m/s]
x0 = [p0;s0]; % state vector

% time parameters
dt = 1; % time step [s]
tf = 10; % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% system matrices
A = [eye(2),eye(2)*dt;zeros(2),eye(2)];
B = [zeros(2);eye(2)*dt];
C = [eye(2),zeros(2)];

% dimension parameters
n = size(A,1);
k = size(C,1);

% noise covariances
Q = eye(2);
R = 9*eye(k);

% simulated noise
w = gaussian_random_sample(zeros(2,1),Q,T-1);
v = gaussian_random_sample(zeros(2,1),R,T-1);

% control signal
u = [-2.5*cos(0.05*t);-2.5*sin(0.05*t)];

% array preallocation
x = zeros(n,T);
y = zeros(k,T);

% stores initial condition
x(:,1) = x0;

% simulation
for tt = 1:(T-1)
    x(:,tt+1) = A*x(:,tt)+B*u(:,tt)+[0;0;w(:,tt)];
    y(:,tt+1) = C*x(:,tt+1)+v(:,tt);
end

% state trajectory with measurement
figure('position',plot_position);
hold on;
plot(x(1,:),x(2,:),'color',cardinal_red,'linewidth',line_width);
plot(y(1,2:end),y(2,2:end),'k:','linewidth',line_width);
hold off;
grid on;
xlabel('$p_{1}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$p_{2}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('true trajectory, $\mathbf{p}$',...
    'position measurement, $\mathbf{y}$','interpreter','latex',...
    'fontsize',legend_font_size);



%% PROBLEM 3

% initial state estimate
mu0 = [1500;100;0;55];
Sigma0 = [250000*eye(2),zeros(2);zeros(2),eye(2)];

% updates process noise covariance to Q'
Q = [zeros(2),zeros(2);zeros(2),eye(2)];

% runs Kalman filter
[mu,Sigma] = KF(A,B,C,Q,R,u,y,mu0,Sigma0);

% plots estimated trajectory with position error ellipses (part (b))
figure('position',plot_position);
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
plot(x(1,:),x(2,:),'k:','linewidth',line_width);
plot(mu(1,2:end),mu(2,2:end),'linewidth',line_width,'color',cardinal_red);
hold off;
grid on;
xlabel('$p_{1}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$p_{2}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('error ellipses for position estimate','true trajectory',...
    'position estimate','interpreter','latex','fontsize',...
    legend_font_size,'location','best');


% plots true trajectory with velocity + velocity error ellipses (part (c))
figure('position',plot_position);
hold on;
for i = 2:length(t)
    [x_ellipse,y_ellipse] = error_ellipse(mu(3:4,i),Sigma(3:4,3:4,i),0.95);
    x_ellipse = x_ellipse+x(1,i);
    y_ellipse = y_ellipse+x(2,i);
    if i == length(t)
        fill(x_ellipse,y_ellipse,'','facecolor',[0.75,0.75,0.75],...
            'edgecolor','none');
        plot(x(1,:),x(2,:),'k:','linewidth',line_width);
        quiver(x(1,i),x(2,i),(10/9)*mu(3,i),(10/9)*mu(4,i),'linewidth',...
            line_width,'color',cardinal_red);
    else
        fill(x_ellipse,y_ellipse,'','facecolor',[0.75,0.75,0.75],...
            'edgecolor','none','handlevisibility','off');
        quiver(x(1,i),x(2,i),(10/9)*mu(3,i),(10/9)*mu(4,i),'linewidth',...
            line_width,'color',cardinal_red,'handlevisibility','off');
    end
end
hold off;
grid on;
xlabel('$p_{1}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$p_{2}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('error ellipses for velocity estimate','true trajectory',...
    'velocity estimate','interpreter','latex','fontsize',...
    legend_font_size,'location','best');



%% PROBLEM 4

% update C matrix
C = [zeros(2),eye(2)];

% new measurement simulation
for tt = 1:(T-1)
    y(:,tt+1) = C*x(:,tt+1)+v(:,tt);
end

% initial state estimate
mu0 = [1000;0;0;50];
Sigma0 = [eye(2),zeros(2);zeros(2),eye(2)];

% runs Kalman filter
[mu,Sigma] = KF(A,B,C,Q,R,u,y,mu0,Sigma0);

% plots estimated trajectory with position error ellipses
figure('position',plot_position);
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
plot(x(1,:),x(2,:),'k:','linewidth',line_width);
plot(mu(1,2:end),mu(2,2:end),'linewidth',line_width,'color',cardinal_red);
hold off;
grid on;
xlabel('$p_{1}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$p_{2}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('error ellipses for position estimate','true trajectory',...
    'position estimate','interpreter','latex','fontsize',...
    legend_font_size,'location','best');



% plots true trajectory with velocity + velocity error ellipses
figure('position',plot_position);
hold on;
for i = 2:length(t)
    [x_ellipse,y_ellipse] = error_ellipse(mu(3:4,i),Sigma(3:4,3:4,i),0.95);
    x_ellipse = x_ellipse+x(1,i);
    y_ellipse = y_ellipse+x(2,i);
    if i == length(t)
        fill(x_ellipse,y_ellipse,'','facecolor',[0.75,0.75,0.75],...
            'edgecolor','none');
        plot(x(1,:),x(2,:),'k:','linewidth',line_width);
        quiver(x(1,i),x(2,i),(10/9)*mu(3,i),(10/9)*mu(4,i),'linewidth',...
            line_width,'color',cardinal_red);
    else
        fill(x_ellipse,y_ellipse,'','facecolor',[0.75,0.75,0.75],...
            'edgecolor','none','handlevisibility','off');
        quiver(x(1,i),x(2,i),(10/9)*mu(3,i),(10/9)*mu(4,i),'linewidth',...
            line_width,'color',cardinal_red,'handlevisibility','off');
    end
end
hold off;
grid on;
xlabel('$p_{1}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
ylabel('$p_{2}\;[\mathrm{m}]$','interpreter','latex','fontsize',...
    axis_font_size);
legend('error ellipses for velocity estimate','true trajectory',...
    'velocity estimate','interpreter','latex','fontsize',...
    legend_font_size,'location','best');