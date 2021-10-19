%==========================================================================
%
% sigma_bounds  Lower and upper sigma bounds for a time series of means and
% covariances.
%
%   [lower_bound,upper_bound] = sigma_bounds(mu,Sigma,N)
%
% Author: Tamas Kis
% Last Update: 2021-08-18
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
function [lower_bound,upper_bound] = sigma_bounds(mu,Sigma,N)
    
    % preallocates matrices
    lower_bound = zeros(size(mu));
    upper_bound = zeros(size(mu));
    
    % finds lower and upper N-sigma bounds
    for i = 1:size(mu,2)
        lower_bound(:,i) = mu(:,i)+N*sqrt(diag(Sigma(:,:,i)));
        upper_bound(:,i) = mu(:,i)-N*sqrt(diag(Sigma(:,:,i)));
    end
    
end