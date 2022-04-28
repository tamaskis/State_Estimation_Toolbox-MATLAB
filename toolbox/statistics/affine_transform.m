%==========================================================================
%
% affine_transform  Affine transformation of mean and covariance.
%
%   [mu_y,Sigma_yy] = affine_transform(mu_x,Sigma_xx,A,b)
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
%   A           - (m×n double) matrix defining linear portion of affine 
%                 transformation Y = AX + b
%   b           - (m×1 double) matrix defining translation portion of
%                 affine transformation Y = AX + b
%
% -------
% OUTPUT:
% -------
%   mu_y        - (m×1 double) mean of Y
%   Sigma_yy    - (m×m double) covariance of Y
%
%==========================================================================
function [mu_y,Sigma_yy] = affine_transform(mu_x,Sigma_xx,A,b)
    
    % mean of Y
    mu_y = A*mu_x+b;
    
    % covariance of Y
    Sigma_yy = A*Sigma_xx*A.';
    
end