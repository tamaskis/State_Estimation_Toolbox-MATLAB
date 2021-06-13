%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Problem Set 2 - Problem 4

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



%% 2D GAUSSIAN #1

% statistics
mu = [1;2];
Sigma = [1,0.5;0.5,2];

% random sample of 1000 points
sample = gaussian_random_sample(mu,Sigma,1000);

% coordinates of error ellipse (P=0.95)
[x,y] = error_ellipse(mu,Sigma,0.95);

% plot
figure('position',plot_position);
hold on;
grid on;
plot(x,y,'linewidth',line_width,'color',cardinal_red);
plot(sample(1,:),sample(2,:),'o','color',[0.5,0.5,0.5]);
xlabel('$x_{1}$','interpreter','latex','fontsize',axis_font_size);
ylabel('$x_{2}$','interpreter','latex','fontsize',axis_font_size);
legend('error ellipse ($P=0.95$)','1000 random samples','interpreter',...
    'latex','fontsize',legend_font_size,'location','best');



%% 2D GAUSSIAN #2

% statistics
mu = [5;2];
Sigma = [1,2;2,7];

% random sample of 1000 points
sample = gaussian_random_sample(mu,Sigma,1000);

% coordinates of error ellipse (P=0.95)
[x,y] = error_ellipse(mu,Sigma,0.95);

% plot
figure('position',plot_position);
hold on;
grid on;
plot(x,y,'linewidth',line_width,'color',cardinal_red);
plot(sample(1,:),sample(2,:),'o','color',[0.5,0.5,0.5]);
xlabel('$x_{1}$','interpreter','latex','fontsize',axis_font_size);
ylabel('$x_{2}$','interpreter','latex','fontsize',axis_font_size);
legend('error ellipse ($P=0.95$)','1000 random samples','interpreter',...
    'latex','fontsize',legend_font_size,'location','best');



%% 2D GAUSSIAN #3

% statistics
mu = [10;10];
Sigma = [2,6;6,20];

% random sample of 1000 points
sample = gaussian_random_sample(mu,Sigma,1000);

% coordinates of error ellipse (P=0.95)
[x,y] = error_ellipse(mu,Sigma,0.95);

% plot
figure('position',plot_position);
hold on;
grid on;
plot(x,y,'linewidth',line_width,'color',cardinal_red);
plot(sample(1,:),sample(2,:),'o','color',[0.5,0.5,0.5]);
xlabel('$x_{1}$','interpreter','latex','fontsize',axis_font_size);
ylabel('$x_{2}$','interpreter','latex','fontsize',axis_font_size);
legend('error ellipse ($P=0.95$)','1000 random samples','interpreter',...
    'latex','fontsize',legend_font_size,'location','best');