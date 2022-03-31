%==========================================================================
%
% EKF  Extended Kalman filter (single iteration).
%
%   [xk,Pk,z_pre,z_post,F_prev,Hk] = EKF(x_prev,P_prev,u_prev,yk,k,fd,...
%       hd,F,H,Q,R)
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
%   yk      - (p×1 double) measurement at current sample time
%   k       - (1×1 double) current sample number
%   fd      - (1×1 function_handle) discrete dynamics equation,
%             xₖ₊₁ = fd(xₖ,uₖ,k) (fd : ℝⁿ×ℝᵐ×ℤ → ℝⁿ)
%   hd      - (1×1 function_handle) discrete measurement equation,
%             yₖ = hd(xₖ,k) (hd : ℝⁿ×ℤ → ℝᵖ)
%   F       - (1×1 function_handle) discrete dynamics Jacobian, 
%             Fₖ = F(xₖ,uₖ,k) (F : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   H       - (1×1 function_handle) discrete measurement Jacobian,
%             Hₖ = H(xₖ,k) (H : ℝⁿ×ℤ → ℝᵖˣⁿ)
%   Q       - (1×1 function_handle) process noise covariance, 
%             Qₖ = Q(xₖ,uₖ,k) (Q : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   R       - (1×1 function_handle) measurement noise covariance, 
%             Rₖ = R(xₖ,k) (R : ℝⁿ×ℤ → ℝᵖˣᵖ)
%
% -------
% OUTPUT:
% -------
%   xk      - (n×1 double) a posteriori state estimate
%   Pk      - (n×n double) a posteriori error covariance
%   z_pre   - (p×1 double) pre-fit measurement residual
%   z_post  - (p×1 double) post-fit measurement residual
%   F_prev  - (n×n double) discrete dynamics Jacobian at previous sample 
%             time
%   Hk      - (p×m double) discrete measurement Jacobian at current sample
%             time
%
%==========================================================================
function [xk,Pk,z_pre,z_post,F_prev,Hk] = EKF(x_prev,P_prev,u_prev,yk,k,...
    fd,hd,F,H,Q,R)
    
    % predict step (time update)
    [x_pred,P_pred,F_prev] = EKF_predict(x_prev,P_prev,u_prev,k,fd,F,Q);

    % update step (measurement update)
    [xk,Pk,z_pre,z_post,Hk] = EKF_update(x_pred,P_pred,yk,k,hd,H,R);
    
end