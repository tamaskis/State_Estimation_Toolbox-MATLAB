%==========================================================================
%
% sigma_bounds  Lower and upper sigma bounds for a time series of means and
% covariances.
%
%   [lower_bound,upper_bound] = sigma_bounds(mu,Sigma,N)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-06-01
%
%--------------------------------------------------------------------------
%
% ------
% INPUTS:
% ------
%   mu      (nxT matrix of nx1 vectors) mean
%   Sigma   (nxnxT array of nxn matrices) covariance
%	N       (1x1) "number of sigmas", i.e. N = 2 returns (+/-)2-sigma
%
% --------
% OUTPUTS:
% --------
%   lower_bound     (nxT matrix of nx1 vectors) lower N-sigma bound
%   upper_bound     (nxT matrix of nx1 vectors) upper N-sigma bound
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