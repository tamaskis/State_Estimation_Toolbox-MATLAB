%% ps2_problem_4
% Problem Set 2, Problem 4
% AA 273 - State Estimation and Filtering for Robotic Perception
%
% Author: Tamas Kis
% Last Update: 2021-08-18



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;

% adds path to root directory and all subdirectories
addpath(genpath("../"));

% loads plot parameters
pp = PLOT_PARAMETERS;



%% PLOTS

% parameters for all plots
N = 1000;   % number of sample points
P = 0.95;   % confidence level

% 2D Gaussian #1
mu = [1;2];
Sigma = [1,0.5;0.5,2];
plot_error_ellipse(mu,Sigma,N,P,pp);

% 2D Gaussian #2
mu = [5;2];
Sigma = [1,2;2,7];
plot_error_ellipse(mu,Sigma,N,P,pp);

% 2D Gaussian #3
mu = [10;10];
Sigma = [2,6;6,20];
plot_error_ellipse(mu,Sigma,N,P,pp);



%% ADDITIONAL FUNCTIONS

%==========================================================================
% plot_error_ellipse  Plot an error ellipse for a 2D Gaussian distribution
% with sample points.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   mu      - (2×1 double) mean
%   Sigma   - (2×2 double) covariance
%   N       - (1×1 double) number of sample points
%   P       - (1×1 double) confidence level (proportion between 0 and 1)
%   pp      - (struct) plot parameters
%
%==========================================================================
function plot_error_ellipse(mu,Sigma,N,P,pp)

    % random sample of N points
    sample = gaussian_random_sample(mu,Sigma,N);

    % coordinates of error ellipse
    [x,y] = error_ellipse(mu,Sigma,P);

    % plot
    figure('position',pp.plot_position);
    hold on;
    grid on;
    plot(x,y,'linewidth',pp.line_width,'color',pp.cardinal_red);
    plot(sample(1,:),sample(2,:),'o','color',[0.5,0.5,0.5]);
    xlabel('$x_{1}$','interpreter','latex','fontsize',pp.axis_font_size);
    ylabel('$x_{2}$','interpreter','latex','fontsize',pp.axis_font_size);
    legend('error ellipse ($P=0.95$)','1000 random samples',...
        'interpreter','latex','fontsize',pp.legend_font_size,'location',...
        'best');

end