%==========================================================================
%
% PLOT_PARAMETERS  Initializes plot parameters.
%
%   plot_parameters = PLOT_PARAMETERS()
%
% Author: Tamas Kis
% Last Update: 2021-08-18
%
%--------------------------------------------------------------------------
%
% -------
% OUTPUT:
% -------
%   pp  - (struct) structure storing plot parameters
%
%==========================================================================
function pp = PLOT_PARAMETERS()

    % plot positions [x,y,l,w]
    plot_position = [540,300,700,500];
    two_subplot_position = [270,300,1200,500];
    three_subplot_position = [0,300,2000,500];
    tall_subplot_position = [540,100,700,800];

    % line width
    line_width = 1.5;
    
    % marker size
    marker_size = 7;

    % font sizes for standard plot
    axis_font_size = 18;
    legend_font_size = 14;
    title_font_size = 18;

    % larger font sizes
    axis_font_size_big = 24;
    legend_font_size_big = 18;
    title_font_size_big = 24;

    % color for plots [rgb]
    cardinal_red = [140,21,21]/255;
    light_cardinal_red = (1-cardinal_red)*0.75+cardinal_red;

    % default MATLAB colors [rgb]
    matlab_blue = [0,0.4470,0.7410];
    matlab_red = [0.8500,0.3250,0.0980];
    matlab_yellow = [0.9290,0.6940,0.1250];
    matlab_purple = [0.4940,0.1840,0.5560];
    matlab_green = [0.4660,0.6740,0.1880];
    matlab_cyan = [0.3010,0.7450,0.9330];
    matlab_maroon = [0.6350,0.0780,0.1840];

    % lighter default MATLAB colors [rgb]
    matlab_light_blue = (1-matlab_blue)*0.75+matlab_blue;
    matlab_light_red = (1-matlab_red)*0.75+matlab_red;
    matlab_light_yellow = (1-matlab_yellow)*0.75+matlab_yellow;
    matlab_light_purple = (1-matlab_purple)*0.75+matlab_purple;
    matlab_light_green = (1-matlab_green)*0.75+matlab_green;
    matlab_light_cyan = (1-matlab_cyan)*0.75+matlab_cyan;
    matlab_light_maroon = (1-matlab_maroon)*0.75+matlab_maroon;

    % packages plot parameters into single structure
    pp.plot_position = plot_position;
    pp.two_subplot_position = two_subplot_position;
    pp.three_subplot_position = three_subplot_position;
    pp.tall_subplot_position = tall_subplot_position;
    pp.line_width = line_width;
    pp.axis_font_size = axis_font_size;
    pp.marker_size = marker_size;
    pp.legend_font_size = legend_font_size;
    pp.title_font_size = title_font_size;
    pp.axis_font_size_big = axis_font_size_big;
    pp.legend_font_size_big = legend_font_size_big;
    pp.title_font_size_big = title_font_size_big;
    pp.cardinal_red = cardinal_red;
    pp.light_cardinal_red = light_cardinal_red;
    pp.matlab_blue = matlab_blue;
    pp.matlab_red = matlab_red;
    pp.matlab_yellow = matlab_yellow;
    pp.matlab_purple = matlab_purple;
    pp.matlab_green = matlab_green;
    pp.matlab_cyan = matlab_cyan;
    pp.matlab_maroon = matlab_maroon;
    pp.matlab_light_blue = matlab_light_blue;
    pp.matlab_light_red = matlab_light_red;
    pp.matlab_light_yellow = matlab_light_yellow;
    pp.matlab_light_purple = matlab_light_purple;
    pp.matlab_light_green = matlab_light_green;
    pp.matlab_light_cyan = matlab_light_cyan;
    pp.matlab_light_maroon = matlab_light_maroon;

end