%==========================================================================
%
% UKF_update  UKF update step (measurement update).
%
%   [xk,Pk,z_pre,z_post] = UKF_update(x_pred,P_pred,yk,k,hd,R)
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
%   k       - (1×1 double) current sample number
%   hd      - (1×1 function_handle) discrete measurement equation,
%             yₖ = hd(xₖ,k) (hd : ℝⁿ×ℤ → ℝᵖ)
%   R       - (1×1 function_handle) measurement noise covariance, 
%             Rₖ = R(xₖ,k) (R : ℝⁿ×ℤ → ℝᵖˣᵖ)
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
function [xk,Pk,z_pre,z_post] = UKF_update(x_pred,P_pred,yk,k,hd,R)
    
    % state (n) and measurement (p) dimensions
    n = length(x_pred);
    p = length(yk);

    % sigma points from predicted state estimate statistics
    [Chi,w] = UT(x_pred,P_pred);
    
    % passing sigma points through nonlinear measurement equation
    Y = zeros(p,2*n+1);
    y_pred = zeros(p,1);
    for i = 1:(2*n+1)
        Y(:,i) = hd(Chi(:,i),k);
        y_pred = y_pred+w(i)*Y(:,i);
    end
    
    % pre-fit measurement residual (innovation)
    z_pre = yk-y_pred;
    
    % covariance of predicted measurement and cross covariance between
    % predicted state and predicted measurement
    Py = zeros(p,p);
    Pxy = zeros(n,p);
    for i = 1:(2*n+1)
        Py = Py+w(i)*(Y(:,i)-y_pred)*(Y(:,i)-y_pred)';
        Pxy = Pxy+w(i)*(Chi(:,i)-x_pred)*(Y(:,i)-y_pred)';
    end
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