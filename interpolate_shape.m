% interpolate_shape  Produces interpolated points along the lines forming a
% shape.
%
%   [x_interp,y_interp] = interpolate_shape(x,y,n)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-24
%
%=========================================================================%
%
% INPUTS:
%   x   x-coordinates of points defining shape
%   y   y-coordinates of points defining shape
%   n   number of equally-spaced interpolated points per line of shape
%
% OUTPUTS:
%   x_new   new x-coordinates of points defining shape
%   y_new   new y-coordinates of points defining shape
%   
%=========================================================================%
function [x_new,y_new] = interpolate_shape(x,y,n)
    
    % preallocates arrays
    x_interp = zeros(n*(length(x)-1),1);
    y_interp = zeros(n*(length(x)-1),1);
    
    for i = 1:(length(x)-1)
        
        % first point forming line
        x1 = x(i);
        y1 = y(i);
        
        % second point forming line
        x2 = x(i+1);
        y2 = y(i+1);
        
        % angle between two points
        theta = atan2(y2-y1,x2-x1);
        
        % distance between two points
        l = sqrt((x2-x1)^2+(y2-y1)^2);
        
        % interpolated points
        for j = 1:n
            x_interp(n*(i-1)+j) = x1+(j/(n+1))*l*cos(theta);
            y_interp(n*(i-1)+j) = y1+(j/(n+1))*l*sin(theta);
        end
        
    end
    
    % "fuses" arrays
    x_new = zeros(length(x)+length(x_interp),1);
    y_new = zeros(size(x_new));
    for i = 1:(length(x)-1)     
        k = 2+(i-1)*(n+1);
        x_new(k-1) = x(i);
        y_new(k-1) = y(i);
        for j = 1:n
            x_new(k+j-1) = x_interp(k+j-i-1);
            y_new(k+j-1) = y_interp(k+j-i-1);
        end
    end
    x_new(end) = x(end);
    y_new(end) = y(end);

end