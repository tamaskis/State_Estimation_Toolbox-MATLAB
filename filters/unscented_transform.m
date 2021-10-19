%==========================================================================
%
% unscented_transform  Unscented transform.
%
%   [Chi,w] = unscented_transform(mu,Sigma,lambda)
%
% Author: Tamas Kis
% Last Update: 2021-08-18
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   mu      - (n×1 double) mean
%   Sigma   - (n×n double) covariance
%   lambda  - (OPTIONAL) (1×1 double) scaling parameter (defaults to 2)
%
% -------
% OUTPUT:
% -------
%   Chi  	- (n×(2n+1) double) matrix of sigma points
%   w       - (1×(2n+1) double) row vector of weights
%
%==========================================================================
function [Chi,w] = unscented_transform(mu,Sigma,lambda)

    % sets default value of lambda to 2
    if (nargin < 3)
        lambda = 2;
    end
    
    % determines dimension of state vector
    n = length(mu);
    
    % preallocates arrays to store sigma points and weights
    Chi = zeros(n,2*n+1);
    w = zeros(1,2*n+1);
    
    % takes the matrix square root used for calculation of the sigma points
    S = chol((lambda+n)*Sigma)';
    
    % obtains the sigma points and their corresponding weights
    Chi(:,1) = mu;
    w(1) = lambda/(lambda+n);
    for i = 2:(n+1)
        Chi(:,i) = mu+S(:,i-1);
        Chi(:,i+n) = mu-S(:,i-1);
        w(i) = 1/(2*(lambda+n));
        w(i+n) = 1/(2*(lambda+n));
    end
    
end