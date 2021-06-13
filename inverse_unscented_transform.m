% inverse_unscented_transform  Inverse unscented transform.
%
%   [mu,Sigma] = inverse_unscented_transform(chi,w)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-22
%
%=========================================================================%
%
% INPUTS:
%   chi     sigma points (n x (2n+1))
%   w       weights (1 x (2n+1))
%
% OUTPUTS:
%   mu      mean (n x 1)
%   Sigma   covariance matrix (n x n)
%
%=========================================================================%
function [mu,Sigma] = inverse_unscented_transform(chi,w)
    
    % determines dimension of state vector
    n = size(chi,1);
    
    % initializes mean and covariance
    mu = zeros(n,1);
    Sigma = zeros(n,n);
    
    % performs the inverse unscented transform
    for i = 1:(2*n+1)
        mu = mu+w(i)*chi(:,i);
    end
    for i = 1:(2*n+1)
        Sigma = Sigma+w(i)*(chi(:,i)-mu)*(chi(:,i)-mu)';
    end
    
end