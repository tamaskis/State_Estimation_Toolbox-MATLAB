%% Copyright (c) 2021 Tamas Kis

% Plot parameters for AA 273 problem sets.



%% PLOT PARAMETERS

% plot positions [x,y,l,w]
plot_position = [540,300,700,500];
two_subplot_position = [270,300,1200,500];
three_subplot_position = [0,300,2000,500];
tall_subplot_position = [540,100,700,800];

% line width [#]
line_width = 1.5;

% font sizes for standard plot [#]
axis_font_size = 18;    % axis label font size
legend_font_size = 14;  % legend font size
title_font_size = 18;   % title font size

% font sizes when plotting for subfigures in LaTeX [#]
axis_font_size_big = 24;    % axis label font size
legend_font_size_big = 18;  % legend font size
title_font_size_big = 24;   % title font size

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