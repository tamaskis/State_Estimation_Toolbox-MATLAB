%==========================================================================
%
% unscented_transform  Unscented transform.
%
%   [chi,W] = unscented_transform(mu,Sigma)
%   [chi,W] = unscented_transform(mu,Sigma,lambda)
%
% Author: Tamas Kis
% Last Update: 2022-04-16
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   mu      - (n×1 double) mean
%   Sigma   - (n×n double) covariance
%   lambda  - (1×1 double) (OPTIONAL) scaling parameter (defaults to 2)
%
% -------
% OUTPUT:
% -------
%   chi     - (n×(2n+1) double) matrix of sigma points
%   W       - ((2n+1)×1 double) vector of weights
%
%==========================================================================
function [chi,W] = unscented_transform(mu,Sigma,lambda)

    % sets default value of lambda to 2
    if (nargin < 3)
        lambda = 2;
    end
    
    % determines dimension of state vector
    n = length(mu);
    
    % preallocates arrays to store sigma points and weights
    chi = zeros(n,2*n+1);
    W = zeros(2*n+1,1);
    
    % square root of covariance matrix (Σ¹ᐟ²) via Cholesky decomposition
    Sigma_sqrt = chol(Sigma)';

    % S matrix
    S = sqrt(lambda+n)*Sigma_sqrt;
    
    % obtains the sigma points and their corresponding weights
    chi(:,1) = mu;
    W(1) = lambda/(lambda+n);
    for i = 2:(n+1)
        chi(:,i) = mu+S(:,i-1);
        chi(:,i+n) = mu-S(:,i-1);
        W(i) = 1/(2*(lambda+n));
        W(i+n) = 1/(2*(lambda+n));
    end
    
end