%==========================================================================
%
% sigma_bounds  Lower and upper sigma bounds for a time series of means and
% covariances.
%
%   [x_lower,x_upper] = sigma_bounds(x,P)
%   [x_lower,x_upper] = sigma_bounds(x,P,M)
%
% Author: Tamas Kis
% Last Update: 2022-04-03
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (n×N double) a posteriori state estimates
%   P       - (n×n×N double) a posteriori error covariances
%   M       - (1×1 double) (OPTIONAL) number of standard deviations for
%             bounds (defaults to 1)
%               --> for example, M = 2 returns ±2σ bounds
%
% -------
% OUTPUT:
% -------
%   x_lower - (n×N double) lower Mσ bound on state estimates
%   x_upper - (n×N double) upper Mσ bound on state estimates
%
% -----
% NOTE:
% -----
%   --> n = state dimension
%   --> N = number of samples
%
%==========================================================================
function [x_lower,x_upper] = sigma_bounds(x,P,M)
    
    % defaults M to 1 if not specified
    if nargin == 2, M = 1; end

    % standard deviations at each sample time
    sigma = covar2stdev(P);

    % finds lower and upper Mσ bounds
    x_lower = x-M*sigma;
    x_upper = x+M*sigma;
    
end