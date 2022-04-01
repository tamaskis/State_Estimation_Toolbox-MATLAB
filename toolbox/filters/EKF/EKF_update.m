%==========================================================================
%
% EKF_update  EKF update step (measurement update).
%
%   [xk,Pk,z_pre,z_post,Hk] = EKF_update(x_pred,P_pred,yk,hd,H,R,k)
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
%   hd      - (1×1 function_handle) discrete measurement equation,
%             yₖ = hd(xₖ,k) (hd : ℝⁿ×ℤ → ℝᵖ)
%   H       - (1×1 function_handle) discrete measurement Jacobian,
%             Hₖ = H(xₖ,k) (H : ℝⁿ×ℤ → ℝᵖˣⁿ)
%   R       - (1×1 function_handle) measurement noise covariance, 
%             Rₖ = R(xₖ,k) (R : ℝⁿ×ℤ → ℝᵖˣᵖ)
%   k       - (1×1 double) current sample number
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
%==========================================================================
function [xk,Pk,z_pre,z_post,Hk] = EKF_update(x_pred,P_pred,yk,hd,H,R,k)
    
    % state dimension
    n = length(x_pred);
    
    % discrete measurement Jacobian at current sample time
    Hk = H(x_pred,k);
    
    % pre-fit measurement residual (innovation)
    z_pre = yk-hd(x_pred,k);
    
    % pre-fit measurement residual covariance (innovation covariance)
    S = Hk*P_pred*Hk.'+R(x_pred,k);
    
    % Kalman gain
    Kk = P_pred*Hk.'/S;

    % a posteriori state estimate at current sample time
    xk = x_pred+Kk*z_pre;
    
    % a posteriori error covariance at current sample time
    Pk = (eye(n)-Kk*Hk)*P_pred;
    
    % post-fit measurement residual
    z_post = yk-hd(xk,k);
    
end