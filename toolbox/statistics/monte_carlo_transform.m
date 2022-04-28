%==========================================================================
%
% monte_carlo_transform  Monte Carlo transformation for passing a Gaussian 
% through a nonlinearity.
%
%   [mu_y,Sigma_yy] = monte_carlo_transform(mu_x,Sigma_xx,f)
%   [mu_y,Sigma_yy] = monte_carlo_transform(mu_x,Sigma_xx,f,N)
%
% Author: Tamas Kis
% Last Update: 2022-04-27
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   mu_x        - (n×1 double) mean of X
%   Sigma_xx    - (n×n double) covariance of X
%   f           - (1×1 function_handle) multivariate, vector-valued 
%                 function, Y = f(X) (f : ℝⁿ → ℝᵐ)
%   N           - (1×1 double) (OPTIONAL) sample size (defaults to 1000)
%
% -------
% OUTPUT:
% -------
%   mu_y        - (m×1 double) mean of Y
%   Sigma_yy    - (m×m double) covariance of Y
%
% -----
% NOTE:
% -----
%   --> This implementation assumes that X is a Gaussian random vector.
%
%==========================================================================
function [mu_y,Sigma_yy] = monte_carlo_transform(mu_x,Sigma_xx,f,N)
    
    % defaults sample size to 1000 if not specified
    if nargin < 4
        N = 1000;
    end
    
    % generate random samples of X
    x = randmvn(mu_x,Sigma_xx,N);
    
    % dimension of Y
    m = length(f(mu_x));
    
    % preallocates matrix to store transformations of X
    y = zeros(m,N);
    
    % transforms X to Y
    for i = 1:N
        y(:,i) = f(x(:,i));
    end
    
    % approximates mean of Y using sample mean
    mu_y = zeros(m,1);
    for i = 1:N
        mu_y = mu_y+y(:,i);
    end
    mu_y = mu_y/N;
    
    % approximates covariance of Y using sample covariance
    Sigma_yy = zeros(m,m);
    for i = 1:N
        Sigma_yy = Sigma_yy+(y(:,i)-mu_y)*(y(:,i)-mu_y).';
    end
    Sigma_yy = Sigma_yy/(N-1);
    
end