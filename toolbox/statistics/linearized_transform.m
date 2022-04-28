%==========================================================================
%
% linearized_transform  Linearized transformation for passing a
% distribution through a nonlinearity.
%
%   [mu_y,Sigma_yy] = linearized_transform(mu_x,Sigma_xx,f)
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
%   J           - (1×1 function_handle) Jacobian of f with respect to x,
%                 J(X) = ∂f/∂X (J : ℝⁿ → ℝᵐˣⁿ)
%
% -------
% OUTPUT:
% -------
%   mu_y        - (m×1 double) mean of Y
%   Sigma_yy    - (m×m double) covariance of Y
%
%==========================================================================
function [mu_y,Sigma_yy] = linearized_transform(mu_x,Sigma_xx,f,J)
    
    % evaluates Jacobian of nonlinear transformation at the mean of x
    A = J(mu_x);
    
    % approximates mean of Y via linearization
    mu_y = f(mu_x);
    
    % approximates covariance of Y via linearization
    Sigma_yy = A*Sigma_xx*A.';
    
end