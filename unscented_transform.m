% unscented_transform  Unscented transform.
%
%   [chi,w] = unscented_transform(mu,Sigma,lambda)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-22
%
%=========================================================================%
%
% INPUTS:
%   mu      mean (n x 1)
%   Sigma   covariance matrix (n x n)
%   lambda  (OPTIONAL) parameter for unscented transform (defaults to 2)
%
% OUTPUTS:
%   chi     sigma points (n x (2n+1))
%   w       weights (1 x (2n+1))
%
%=========================================================================%
function [chi,w] = unscented_transform(mu,Sigma,lambda)

    % sets default value of lambda to 2
    if (nargin < 3) || isempty(lambda)
        lambda = 2;
    end
    
    % determines dimension of state vector
    n = length(mu);
    
    % preallocates arrays to store sigma points and weights
    chi = zeros(n,2*n+1);
    w = zeros(1,2*n+1);
    
    % takes the matrix square root used for calculation of the sigma points
    matrix_square_root = chol((lambda+n)*Sigma)';
    
    % performs the unscented transform
    chi(:,1) = mu;
    w(1) = lambda/(lambda+n);
    for i = 2:(n+1)
        chi(:,i) = mu+matrix_square_root(:,i-1);
        chi(:,i+n) = mu-matrix_square_root(:,i-1);
        w(i) = 1/(2*(lambda+n));
        w(i+n) = 1/(2*(lambda+n));
    end
    
end