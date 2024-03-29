%==========================================================================
%
% UKF_update  UKF update step (measurement update).
%
%   [xk,Pk,z_pre,z_post] = UKF_update(x_pred,P_pred,yk,hd,R,k)
%
% Author: Tamas Kis
% Last Update: 2022-04-16
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
%
%==========================================================================
function [xk,Pk,z_pre,z_post] = UKF_update(x_pred,P_pred,yk,hd,R,k)
    
    % function handle for nonlinearity
    g = @(x) hd(x,k);
    
    % predicted measurement, uncorrected measurement error covariance, and 
    % state/measurement cross covariance at current sample time
    [y_pred,Py,Pxy] = unscented_transform(x_pred,P_pred,g,true);
    
    % pre-fit measurement residual (innovation)
    z_pre = yk-y_pred;
    
    % adds effect of measurement noise to correct the measurement error 
    % covariance
    Py = Py+R(x_pred,k);
    
    % Kalman gain
    K = Pxy/Py;
    
    % a posteriori state estimate at current sample time
    xk = x_pred+K*z_pre;
    
    % a posteriori error covariance at current sample time
    Pk = P_pred-K*Py*K.';
    
    % post-fit measurement residual
    z_post = yk-hd(xk,k);
    
end