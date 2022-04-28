%==========================================================================
%
% scaled_unscented_transform  Scaled unscented transformation for passing a
% distribution through a nonlinearity.
%
%   [mu_y,Sigma_yy] = scaled_unscented_transform(mu_x,Sigma_xx,f)
%   [mu_y,Sigma_yy,Sigma_xy] = scaled_unscented_transform(mu_x,Sigma_xx,...
%       f,true)
%   [__] = scaled_unscented_transform(__,alpha,beta,kappa)
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
%   alpha       - (1×1 double) (OPTIONAL) spread parameter, α (defaults to
%                 10⁻³)
%   beta        - (1×1 double) (OPTIONAL) distribution parameter, β 
%                 (defaults to 2, assuming X is Gaussian)
%   kappa       - (1×1 double) (OPTIONAL) secondary scaling parameter, κ 
%                 (defaults to 0)
%
% -------
% OUTPUT:
% -------
%   mu_y        - (m×1 double) mean of Y
%   Sigma_yy    - (m×m double) covariance of Y
%   Sigma_xy    - (n×m double) cross covariance of X and Y
%
%==========================================================================
function [mu_y,Sigma_yy,Sigma_xy] = scaled_unscented_transform(mu_x,...
    Sigma_xx,f,cross_covar,alpha,beta,kappa)
    
    % defaults "alpha" to 10⁻³ if not input
    if (nargin < 5) || isempty(alpha)
        alpha = 1e-3;
    end

    % defaults "beta" to 2 if not input
    if (nargin < 6) || isempty(beta)
        beta = 2;
    end

    % defaults "kappa" to 0 if not input
    if (nargin < 7) || isempty(kappa)
        kappa = 0;
    end
    
    % defaults "cross_covar" to "false" if not input
    if (nargin < 4) || isempty(cross_covar)
        cross_covar = false;
    end
    
    % ------------------------
    % Generating sigma points.
    % ------------------------
    
    % dimension of X
    n = length(mu_x);
    
    % square root of covariance matrix (Σ¹ᐟ²) via Cholesky decomposition
    Sigma_sqrt = chol(Sigma_xx)';
    
    % scaling parameter
    lambda = alpha^2*(n+kappa)-n;
    
    % S matrix
    S = sqrt(n+lambda)*Sigma_sqrt;
    
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
    w0 = lambda/(n+lambda);
    w = 1/(2*(n+lambda));
    
    % γ parameter
    gamma = 1-alpha^2+beta;
    
    % approximates mean of Y using sample mean
    mu_y =  w0*y0;
    for i = 1:(2*n)
        mu_y = mu_y+w*y(:,i);
    end
    
    % approximates covariance of Y using sample covariance
    Sigma_yy = (w0+gamma)*(y0-mu_y)*(y0-mu_y).';
    for i = 1:(2*n)
        Sigma_yy = Sigma_yy+w*(y(:,i)-mu_y)*(y(:,i)-mu_y).';
    end
    
    % approximates cross covariance of X and Y using sample cross covar.
    if cross_covar
        Sigma_xy = zeros(n,m);
        for i = 1:(2*n)
            Sigma_xy = Sigma_xy+w*(x(:,i)-mu_x)*(y(:,i)-mu_y).';
        end
    else
        Sigma_xy = NaN;
    end
    
end