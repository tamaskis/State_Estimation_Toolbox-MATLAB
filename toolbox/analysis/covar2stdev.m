%==========================================================================
%
% covar2stdev  Standard deviations from covariance matrix.
%
%   sigma = covar2stdev(Sigma)
%
% Author: Tamas Kis
% Last Update: 2022-04-03
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   Sigma   - (n×n×N double) time history of covariances (Σ)
%
% -------
% OUTPUT:
% -------
%   sigma   - (n×N double) time history of standard deviations (σ)
%
% -----
% NOTE:
% -----
%   --> n = random vector dimension
%   --> N = number of samples
%
%==========================================================================
function sigma = covar2stdev(Sigma)
    
    % random vector dimension (n) and number of samples (N)
    n = size(Sigma,1);
    N = size(Sigma,3);

    % preallocates array to store n standard deviations at N sample times
    sigma = zeros(n,N);
    
    % calculates standard deviation of each random variable at each sample
    for i = 1:n
        for k = 1:N
            sigma(i,k) = sqrt(Sigma(i,i,k));
        end
    end
    
end