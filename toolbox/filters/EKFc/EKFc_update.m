%==========================================================================
%
% EKFc_update  EKF update step (measurement update) for continuous-time
% systems.
%
%   [xk,Pk,z_pre,z_post,Hk] = EKFc_update(x_pred,P_pred,yk,h,C,R,k,dt,t0)
%
% Author: Tamas Kis
% Last Update: 2022-03-31
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x_pred  - (n×1 double) a priori state estimate at current sample time
%   P_pred  - (n×n double) a priori error covariance at current sample
%   yk      - (p×1 double) measurement at current sample time
%   h       - (1×1 function_handle) continuous measurement equation,
%             y = h(x,t) (h : ℝⁿ×ℝ → ℝᵖ)
%   C       - (1×1 function_handle) continuous measurement Jacobian,
%             C(t) = C(x,t) (C : ℝⁿ×ℝ → ℝᵖˣⁿ)
%   R       - (1×1 function_handle) measurement noise covariance, 
%             Rₖ = R(xₖ,k) (R : ℝⁿ×ℤ → ℝᵖˣᵖ)
%   k       - (1×1 double) current sample number
%   dt      - (1×1 double) time step, Δt
%   t0      - (1×1 double) initial time, t₀
%
% -------
% OUTPUT:
% -------
%   xk      - (n×1 double) a posteriori state estimate at current sample 
%             time
%   Pk      - (n×n double) a posteriori error covariance at current sample
%             time
%   z_pre   - (p×1 double) pre-fit measurement residual
%   z_post  - (p×1 double) post-fit measurement residual
%   Hk      - (p×m double) discrete measurement Jacobian at current sample
%             time
%
% -----
% NOTE:
% -----
%   --> This function makes use of the fact that Hₖ = C(tₖ).
%
%==========================================================================
function [xk,Pk,z_pre,z_post,Hk] = EKFc_update(x_pred,P_pred,yk,h,C,R,k,...
    dt,t0)
    
    % current sample time
    tk = k2t_num(k,dt,t0);

    % state dimension
    n = length(x_pred);
    
    % discrete measurement Jacobian at current sample time
    Hk = C(x_pred,tk);

    % measurement noise covariance at current sample time
    Rk = R(x_pred,k);
    
    % pre-fit measurement residual (innovation)
    z_pre = yk-h(x_pred,tk);
    
    % pre-fit measurement residual covariance (innovation covariance)
    S = Hk*P_pred*Hk.'+R(x_pred,k);
    
    % Kalman gain
    Kk = P_pred*Hk.'/S;

    % auxiliary matrix for use with Joseph's formula
    Jk = eye(n)-Kk*Hk;
    
    % a posteriori state estimate at current sample time
    xk = x_pred+Kk*z_pre;
    
    % a posteriori error covariance at current sample time (Joseph's
    % formula)
    Pk = Jk*P_pred*Jk.'+Kk*Rk*Kk.';
    
    % post-fit measurement residual
    z_post = yk-h(xk,tk);
    
end