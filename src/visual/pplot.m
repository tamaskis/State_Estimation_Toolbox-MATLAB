%==========================================================================
%
% pplot  TODO.
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
%   pp          - (struct) TODO
%   x           - (n×1 or 1×n double) independent variable data
%   y           - (n×1 or 1×n double) dependent variable data
%   settings    - (OPTIONAL) (struct) TODO
%
%==========================================================================
function pplot(pp,x,y,settings)

    % logicals that dictate which settings need to be defaulted
    default_line_color = ~((nargin == 4) && (isfield(settings,'color')));
    default_line_width = ~((nargin == 4) && (isfield(settings,...
        'line_width')));
    default_grid = ~((nargin == 4) && (isfield(settings,'grid')));
    default_position = ~((nargin == 4) && (isfield(settings,'position')));
    
    % initializes figure
    if default_position
        figure('position',pp.two_subplot_position);
    else
    end
    
    % sets line color
    
    
    % sets line width
    if default_line_width
        line_width = pp.line_width;
    else
        line_width = settings.line_width;
    end
    
    % initializes figure
    figure('position',pp.two_subplot_position);

    % plots line
    plot(x,y,'color',pp.cardinal_red,'linewidth',line_width);
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