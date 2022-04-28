%==========================================================================
%
% linear_transform  Linear transformation of mean and covariance.
%
%   [mu_y,Sigma_yy] = linear_transform(mu_x,Sigma_xx,A)
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
%   A           - (m×n double) matrix defining linear transformation Y = AX
%
% -------
% OUTPUT:
% -------
%   mu_y        - (m×1 double) mean of Y
%   Sigma_yy    - (m×m double) covariance of Y
%
%==========================================================================
function [mu_y,Sigma_yy] = linear_transform(mu_x,Sigma_xx,A)
    
    % mean of Y
    mu_y = A*mu_x;
    
    % covariance of Y
    Sigma_yy = A*Sigma_xx*A.';
    
end