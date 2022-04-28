%==========================================================================
%
% unscented_transform  Unscented transformation for passing a Gaussian 
% through a nonlinearity.
%
%   [mu_y,Sigma_yy] = unscented_transform(mu_x,Sigma_xx,f)
%   [mu_y,Sigma_yy,Sigma_xy] = unscented_transform(mu_x,Sigma_xx,f,true)
%   [__] = unscented_transform(__,lambda)
%
% Author: Tamas Kis
% Last Update: 2022-04-16
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
%   cross_covar - (1×1 logical) (OPTIONAL) input as "true" if cross 
%                 covariance of X and Y should be calculated) (defaults to 
%                 false)
%   lambda      - (1×1 double) (OPTIONAL) scaling parameter, λ (defaults to
%                 λ = 2)
%
% -------
% OUTPUT:
% -------
%   mu_y        - (m×1 double) mean of Y
%   Sigma_yy    - (m×m double) covariance of Y
%   Sigma_xy    - (n×m double) cross covariance of X and Y
%
%==========================================================================
function [mu_y,Sigma_yy,Sigma_xy] = unscented_transform(mu_x,Sigma_xx,...
    cross_covar,lambda)
    
    % defaults "cross_covar" to "false" if not input
    if (nargin < 3) || isempty(cross_covar)
        cross_covar = false;
    end
    
    % defaults "lambda" to 2 if not input
    if (nargin < 4) || isempty(lambda)
        lambda = 2;
    end
    
    % ------------------------
    % Generating sigma points.
    % ------------------------
    
    % dimension of X
    n = length(mu_x);
    
    % preallocates arrays to store sigma points and weights
    x = zeros(n,2*n+1);
    W = zeros(2*n+1,1);
    
    % square root of covariance matrix (Σ¹ᐟ²) via Cholesky decomposition
    Sigma_sqrt = chol(Sigma_xx)';
    
    % S matrix
    S = sqrt(lambda+n)*Sigma_sqrt;
    
    % obtains the sigma points and their corresponding weights
    x(:,1) = mu_x;
    W(1) = lambda/(lambda+n);
    for i = 2:(n+1)
        x(:,i) = mu_x+S(:,i-1);
        x(:,i+n) = mu_x-S(:,i-1);
        W(i) = 1/(2*(lambda+n));
        W(i+n) = 1/(2*(lambda+n));
    end
    
    % ------------------------------------------
    % Passing sigma points through nonlinearity.
    % ------------------------------------------
    
    % passes 1st sigma point through nonlinearity
    y1 = f(x(:,1));
    
    % dimension of Y
    m = length(Y1);
    
    % preallocates array to store samples of Y
    y = zeros(m,2*n+1);
    
    % stores first sample of Y
    y(:,1) = y1;
    
    % passes sigma points through the nonlinearity
    for i = 2:(2*n+1)
        y(:,i) = f(x(:,i));
    end

    % -------------------------------------
    % Calculating mean and covariance of Y.
    % -------------------------------------
    
    % approximates mean of Y using sample mean
    mu_y = zeros(m,1);
    for i = 1:(2*n+1)
        mu_y = mu_y+W(i)*y(:,i);
    end
    
    % approximates covariance of Y using sample covariance
    Sigma_yy = zeros(m,m);
    for i = 1:(2*n+1)
        Sigma_yy = W(i)*Sigma_yy+(y(:,i)-mu_y)*(y(:,i)-mu_y).';
    end
    
    % ----------------------------------------
    % Calculating cross covariance of X and Y.
    % ----------------------------------------

    % approximates using sample cross covariance
    if cross_covar
        Sigma_xy = zeros(n,m);
        for i = 1:(2*n+1)
            Sigma_xy = Sigma_xy+W(i)*(x(:,i)-mu_x)*(y(:,i)-mu_y).';
        end
    end
    
end