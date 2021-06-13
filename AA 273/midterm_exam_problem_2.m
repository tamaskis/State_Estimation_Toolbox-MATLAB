%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Midterm Exam Problem 2

%-------------------------------------------------------------------------%



%% SCRIPT SETUP

% clears variables and command window, closes all figures
clear;
clc;
close all;

% plot parameters
plot_position = [540,300,700,500]; % plot position [x,y,l,w]
line_width = 1.5; % line width [#]
axis_font_size = 18; % axis label font size [#]
legend_font_size = 14; % legend font size [#]
title_font_size = 20; % title font size [#]



%% PART (E)

% initial prior
mu0 = [0;0];
Sigma0 = eye(2);

% noise variance
sigma_n_2 = 2;

% training data
x_data = [0,-2,1,-1,2];
y_data = [0.438,4.28,-0.811,2.70,-2.39];

% arrays to store posterior distributions
mu = zeros(2,6);
Sigma = zeros(2,2,6);

% assigns initial conditions to arrayss
mu(:,1) = mu0;
Sigma(:,:,1) = Sigma0;

% recursive estimation
for i = 2:6
    
    % gets x and y from data
    x = [1;x_data(i-1)];
    y = y_data(i-1);
    
    % M and b
    M = inv(Sigma(:,:,i-1))+(x*x')/sigma_n_2;
    b = (y*x/sigma_n_2)+inv(Sigma(:,:,i-1))*mu(:,i-1);
    
    % posterior distribution (updates mean and covariance)
    mu(:,i) = inv(M)*b;
    Sigma(:,:,i) = inv(M);
    
end

% calculates posterior predictive distribution
xstar = [ones(size(-2:0.1:2));-2:0.1:2];
imax = size(xstar,2);
mean_predict = zeros(6,imax);
for i = 1:imax
    for j = 1:6
        mean_predict(j,i) = xstar(:,i)'*mu(:,j);
    end
end

% posterior predictive distribution mean
figure('position',plot_position);
hold on;
for i = 1:6
    plot(xstar(2,:),mean_predict(i,:),'o','markersize',7,'linewidth',...
        line_width);
end
hold off;
grid on;
xlabel('$x_{*}$','interpreter','latex','fontsize',axis_font_size);
ylabel("$\mu_{y_{*}}=\mathbf{x}_{*}^{T}$\boldmath$\mu$\unboldmath$'$",...
    'interpreter','latex','fontsize',axis_font_size);
legend('prior','$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','interpreter',...
    'latex','fontsize',legend_font_size,'location','best');

% posterior distribution error ellipses
figure('position',plot_position);
hold on;
for i = 1:6
    [x_ellipse,y_ellipse] = error_ellipse(mu(:,i),Sigma(:,:,i),0.95);
    plot(x_ellipse,y_ellipse,'linewidth',line_width);
end
hold off;
grid on;
axis equal;
xlabel('$\theta_{1}$','interpreter','latex','fontsize',axis_font_size);
ylabel('$\theta_{2}$','interpreter','latex','fontsize',axis_font_size);
legend('prior','$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','interpreter',...
    'latex','fontsize',legend_font_size,'location','best');



%% ADDITIONAL FUNCTIONS

% INPUT: mu - mean (2-by-1 vector)
%        Sigma - covariance matrix (2-by-2 matrix)
%        P - confidence level as proportion (scalar between 0 and 1)
% OUTPUT: [x,y] - coordinates of the error ellipse (629-by-1 vectors)
function [x,y] = error_ellipse(mu,Sigma,P)
    
    % domain for theta
    theta  = 0:0.01:2*pi;
    
    % radius of circle
    epsilon = (1-P)/(2*pi*sqrt(det(Sigma)));
    r0 = sqrt(-2*log(2*pi*epsilon*sqrt(det(Sigma))));
    
    % defines circle of radius r0 located at origin
    We = r0*[cos(theta);sin(theta)];
    
    % transforms w to x
    Xe = sqrtm(Sigma)*We+mu;

    % extracts coordinates of error ellipse
    x = Xe(1,:)';
    y = Xe(2,:)';
    
end