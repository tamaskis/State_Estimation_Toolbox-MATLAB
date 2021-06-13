%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Problem Set 2 - Problem 5

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



%% SOLUTION

% measurement
y = [0.27;0.62];

% statistics
mu_x = [0;0];
Sigma_xx = [0.7,0.73;0.73,1.1];
mu_y = [0;0];
Sigma_yy = [0.7,0.19;0.19,0.16];
Sigma_xy = [0.63,0.23;0.72,0.31];

% posterior mean and covariance of rover position
mu_xgy = mu_x+Sigma_xy*inv(Sigma_yy)*(y-mu_y)
Sigma_xgy = Sigma_xx-Sigma_xy*inv(Sigma_yy)*Sigma_xy'

% coordinates of error ellipses (P=0.95)
[x_prior,y_prior] = error_ellipse(mu_x,Sigma_xx,0.95);
[x_posterior,y_posterior] = error_ellipse(mu_xgy,Sigma_xgy,0.95);

% plot
figure('position',plot_position);
hold on;
grid on;
plot(x_prior,y_prior,'linewidth',line_width);
plot(x_posterior,y_posterior,'linewidth',line_width);
xlabel('$x_{1}$','interpreter','latex','fontsize',axis_font_size);
ylabel('$x_{2}$','interpreter','latex','fontsize',axis_font_size);
legend('prior','posterior','interpreter','latex','fontsize',...
    legend_font_size,'location','best');