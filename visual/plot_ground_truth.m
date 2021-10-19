%==========================================================================
%
% plot_ground_truth  Produces plots of the ground truth and its
% measurement.
%
%   TODO
%
% Author: Tamas Kis
% Last Update: 2021-08-23
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   mu              - (n×T double) time history of mean
%   Sigma           - (n×n×T double) time history of covariance
%	N               - (1×1 double) "number of sigmas"
%                       --> for example, N = 2 returns (+/-)2-sigma
%
% -------
% OUTPUT:
% -------
%   lower_bound     - (n×T double) lower N-sigma bound
%   upper_bound     - (n×T double) upper N-sigma bound
%
% -----
% NOTE:
% -----
%   --> "mu", "lower_bound", and "upper_bound" are n×T matrices of n×1
%       vectors.
%   --> "Sigma" is an n×n×T array of n×n matrices.
%
%==========================================================================
function plot_ground_truth(pp,x,y,varargin)
    
    % initializes figure
    figure('position',pp.two_subplot_position);

    % true state
    subplot(1,2,1);
    plot(t,x(1,:),'color',pp.cardinal_red,'linewidth',pp.line_width);
    grid on;
    axis square;
    xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
        pp.axis_font_size);
    ylabel('$\omega_{x}\;[\mathrm{rad/s}]$','interpreter','latex',...
        'fontsize', pp.axis_font_size);
    title('\textbf{True State}','interpreter','latex','fontsize',...
        pp.title_font_size)

    % measurement
    subplot(1,2,2);
    plot(t(2:end),y(1,2:end),'color',pp.cardinal_red,'linewidth',...
        pp.line_width);
    grid on;
    axis square;
    xlabel('$t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
        pp.axis_font_size);
    ylabel('$\omega_{x}\;[\mathrm{rad/s}]$','interpreter','latex',...
        'fontsize', pp.axis_font_size);
    title('\textbf{Measurement}','interpreter','latex','fontsize',...
        pp.title_font_size)
    
end