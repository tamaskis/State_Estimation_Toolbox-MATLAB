% error_ellipse  Draws an error ellipse for a two-dimensional Gaussian
% random vector.
%
%   [x,y] = error_ellipse(mu,Sigma,P)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-22
%
%=========================================================================%
%
% INPUTS:
%   mu      mean (2 x 1)
%   Sigma   covariance (2 x 2)
%	P       confidence level as a proportion (scalar between 0 and 1)
%
% OUTPUTS:
%   x       x-coordinates of ellipse (629 x 1)
%   y       y-coordinates of ellipse (629 x 1)
%
%=========================================================================%
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