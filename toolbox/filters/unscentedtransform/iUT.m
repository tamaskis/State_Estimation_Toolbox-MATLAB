%==========================================================================
%
% iUT  Inverse unscented transform.
%
%   [mu,Sigma] = iUT(Chi,w)
%
% Author: Tamas Kis
% Last Update: 2021-11-18
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   Chi     - (n×(2n+1) double) matrix of sigma points
%   w       - (1×(2n+1) double) row vector of weights
%
% -------
% OUTPUT:
% -------
%   mu      - (n×1 double) mean
%   Sigma   - (n×n double) covariance
%
%==========================================================================
function [mu,Sigma] = iUT(Chi,w)
    
    % determines dimension of state vector
    n = size(Chi,1);
    
    % initializes mean and covariance
    mu = zeros(n,1);
    Sigma = zeros(n,n);
    
    % approximates mean with the sample mean
    for i = 1:(2*n+1)
        mu = mu+w(i)*Chi(:,i);
    end
    
    % approximates covariance with the sample covariance
    for i = 1:(2*n+1)
        Sigma = Sigma+w(i)*(Chi(:,i)-mu)*(Chi(:,i)-mu)';
    end
    
end