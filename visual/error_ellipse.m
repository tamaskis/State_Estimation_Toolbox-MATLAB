%==========================================================================
%
% error_ellipse  Draws an error ellipse for a two-dimensional Gaussian
% random vector.
%
%   [x,y] = error_ellipse(mu,Sigma,P)
%
% Last Update: 2021-7-26
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   mu      - (2×1) mean
%   Sigma   - (2×2) covariance
%	P       - (1×1) confidence level as a proportion (between 0 and 1)
%
% -------
% OUTPUT:
% -------
%   x       - (629×1) x-coordinates of ellipse
%   y       - (629×1) y-coordinates of ellipse
%
%==========================================================================
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