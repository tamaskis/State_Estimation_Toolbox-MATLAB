%==========================================================================
%
% unscented_transform  Unscented transformation for passing a distribution
% through a nonlinearity.
%
%   [mu_y,Sigma_yy] = unscented_transform(mu_x,Sigma_xx,f)
%   [mu_y,Sigma_yy,Sigma_xy] = unscented_transform(mu_x,Sigma_xx,f,true)
%   [__] = unscented_transform(__,kappa)
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
%   cross_covar - (1×1 logical) (OPTIONAL) input as "true" if cross 
%                 covariance of X and Y should be calculated) (defaults to 
%                 false)
%   kappa       - (1×1 double) (OPTIONAL) scaling parameter, κ (defaults to
%                 κ = 3-n)
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
    cross_covar,kappa)
    
    % dimension of X
    n = length(mu_x);
    
    % defaults "kappa" to 3-n if not input
    if (nargin < 4) || isempty(kappa)
        kappa = 3-n;
    end
    
    % defaults "cross_covar" to "false" if not input
    if (nargin < 3) || isempty(cross_covar)
        cross_covar = false;
    end
    
    % ------------------------
    % Generating sigma points.
    % ------------------------
    
    % square root of covariance matrix (Σ¹ᐟ²) via Cholesky decomposition
    Sigma_sqrt = chol(Sigma_xx)';
    
    % S matrix
    S = sqrt(n+kappa)*Sigma_sqrt;
    
    % n×n matrix storing mean of X in each column
    M = repmat(mu_x,1,n);
    
    % sigma point matrix
    x = [(M+S),(M-S)];
    
    % ------------------------------------------
    % Passing sigma points through nonlinearity.
    % ------------------------------------------
    
    % passes 1st sigma point through nonlinearity
    y0 = f(mu_x);
    
    % dimension of Y
    m = length(y0);
    
    % preallocates array to store samples of Y
    y = zeros(m,2*n);
    
    % passes remaining sigma points through the nonlinearity
    for i = 1:(2*n)
        y(:,i) = f(x(:,i));
    end
    
    % ----------------------------------------------------------
    % Mean and covariance of Y, and cross covariance of X and Y.
    % ----------------------------------------------------------
    
    % determines the weights
    w0 = kappa/(n+kappa);
    w = 1/(2*(n+kappa));
    
    % approximates mean of Y using sample mean
    mu_y =  w0*y0;
    for i = 1:(2*n)
        mu_y = mu_y+w*y(:,i);
    end
    
    % approximates covariance of Y using sample covariance
    if (kappa >= 0)
        Sigma_yy = w0*(y0-mu_y)*(y0-mu_y).';
    else
        Sigma_yy = zeros(m,m);
    end
    for i = 1:(2*n)
        Sigma_yy = Sigma_yy+w*(y(:,i)-mu_y)*(y(:,i)-mu_y).';
    end
    
    % approximates cross covariance of X and Y using sample cross covar.
    if cross_covar
        Sigma_xy = zeros(n,m);
        for i = 1:(2*n)
            Sigma_xy = Sigma_xy+w*(x(:,i)-mu_x)*(y(:,i)-mu_y).';
        end
    end
    
end