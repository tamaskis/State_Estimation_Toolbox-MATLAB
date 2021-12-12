% gaussian_random_sample  Returns a sample of N vectors of (stored as an
% n-by-N matrix) of an n-dimensional Gaussian random vector X.
%
%   X = gaussian_random_sample(mu,Sigma,N)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-22
%
%=========================================================================%
%
% INPUTS:
%   mu      mean (n x 1)
%   Sigma   covariance (n x n)
%	N       (OPTIONAL) number of samples (defaults to 1)
%
% OUTPUTS:
%   X       samples of random vector, each column stores one realization/
%           sample (n x N)
%
%=========================================================================%
function X = gaussian_random_sample(mu,Sigma,N)

    % sets default number of samples to 1
    if (nargin < 3) || isempty(N)
        N = 1;
    end
    
    % determines dimension of random vector
    n = length(mu);
    
    % N random samples from standard normal distribution
    W = randn(n,N);
    
    % square root of covariance matrix using Cholesky decomposition
    Sigma_sqrt = chol(Sigma)';
    
    % n-by-N matrix M storing mean in each column
    M = repmat(mu,1,N);

    % converts to Gaussian distribution with given parameters
    X = Sigma_sqrt*W+M;
    
end