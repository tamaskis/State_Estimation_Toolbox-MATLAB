%==========================================================================
%
% UKF_predict  UKF predict step (time update).
%
%   [x_pred,P_pred] = UKF_predict(x_prev,P_prev,u_prev,fd,Q,k)
%
% Author: Tamas Kis
% Last Update: 2022-03-31
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x_prev  - (n×1 double) state estimate at previous sample time
%   P_prev  - (n×n double) error covariance at previous sample time
%   u_prev  - (m×1 double) control input at previous sample time
%   fd      - (1×1 function_handle) discrete dynamics equation,
%             xₖ₊₁ = fd(xₖ,uₖ,k) (fd : ℝⁿ×ℝᵐ×ℤ → ℝⁿ)
%   Q       - (1×1 function_handle) process noise covariance, 
%             Qₖ = Q(xₖ,uₖ,k) (Q : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   k       - (1×1 double) current sample number
%
% -------
% OUTPUT:
% -------
%   x_pred  - (n×1 double) a priori state estimate at current sample time
%   P_pred  - (n×n double) a priori error covariance at current sample time
%
%==========================================================================
function [x_pred,P_pred] = UKF_predict(x_prev,P_prev,u_prev,fd,Q,k)

    % function handle for nonlinearity
    g = @(x) fd(x,u_prev,k-1);

    % a priori state estimate and uncorrected error covariance at current 
    % sample time
    [x_pred,P_tilde] = unscented_transform(x_prev,P_prev,g);
    
    % adds effect of process noise to correct the error covariance
    P_pred = P_tilde+Q(x_prev,u_prev,k-1);
    
end