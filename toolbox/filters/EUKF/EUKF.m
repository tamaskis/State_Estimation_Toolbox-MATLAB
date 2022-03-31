%==========================================================================
%
% EUKF  Extended/Unscented Kalman filter.
%
%   [xk,Pk,z_pre,z_post] = EUKF(x_prev,P_prev,u_prev,yk,k,fd,hd,F,Q,R)
%
% Author: Tamas Kis
% Last Update: 2022-03-28
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x_prev  - (n×1 double) state estimate at previous sample time
%   P_prev  - (n×n double) error covariance at previous sample time
%   u_prev  - (m×1 double) control input at previous sample time
%   yk      - (p×1 double) measurement at current sample time
%   k       - (1×1 double) current sample number
%   fd      - (1×1 function_handle) discrete nonlinear dynamics equation,
%             xₖ₊₁ = fd(xₖ,uₖ,k) (fd : ℝⁿ×ℝᵐ×ℤ → ℝⁿ)
%   hd      - (1×1 function_handle) discrete nonlinear measurement 
%             equation, yₖ = hd(xₖ,k) (fd : ℝⁿ×ℤ → ℝᵖ)
%   F       - (1×1 function_handle) Fₖ = F(xₖ,uₖ,k) --> discrete dynamics
%             Jacobian (F : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   Q       - (1×1 function_handle) Qₖ = Q(xₖ,uₖ,k) --> process noise 
%             covariance (Q : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   R       - (1×1 function_handle) Rₖ = R(xₖ,k) --> measurement noise 
%             covariance (R : ℝⁿ×ℤ → ℝᵖˣᵖ)
%
% -------
% OUTPUT:
% -------
%   xk      - (n×1 double) a posteriori state estimate
%   Pk      - (n×n double) a posteriori error covariance
%   z_pre   - (p×1 double) pre-fit measurement residual
%   z_post  - (p×1 double) post-fit measurement residual
%
%==========================================================================
function [xk,Pk,z_pre,z_post] = EUKF(x_prev,P_prev,u_prev,yk,k,fd,hd,F,Q,R)
    
    % predict step (time update)
    [x_pred,P_pred] = EKF_predict(x_prev,P_prev,u_prev,k,fd,F,Q);
    
    % update step (measurement update)
    [xk,Pk,z_pre,z_post] = UKF_update(x_pred,P_pred,yk,k,hd,R);
    
end