%==========================================================================
%
% EKF_predict  EKF predict step (time update).
%
%   [x_pred,P_pred,F_prev] = EKF_predict(x_prev,P_prev,u_prev,fd,F,Q,k)
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
%   F       - (1×1 function_handle) discrete dynamics Jacobian, 
%             Fₖ = F(xₖ,uₖ,k) (F : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   Q       - (1×1 function_handle) process noise covariance, 
%             Qₖ = Q(xₖ,uₖ,k) (Q : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   k       - (1×1 double) current sample number
%
% -------
% OUTPUT:
% -------
%   x_pred  - (n×1 double) a priori state estimate at current sample time
%   P_pred  - (n×n double) a priori error covariance at current sample time
%   F_prev  - (n×n double) discrete dynamics Jacobian at previous sample 
%             time
%
%==========================================================================
function [x_pred,P_pred,F_prev] = EKF_predict(x_prev,P_prev,u_prev,fd,F,...
    Q,k)
    
    % discrete dynamics Jacobian at previous sample time
    F_prev = F(x_prev,u_prev,k-1);
    
    % a priori state estimate at current sample time
    x_pred = fd(x_prev,u_prev,k-1);
    
    % a priori error covariance at current sample time
    P_pred = F_prev*P_prev*F_prev.'+Q(x_prev,u_prev,k-1);
    
end