%==========================================================================
%
% EKF  Extended Kalman filter.
%
%   [x,P,tsol,rank_Ob,z_pre,z_post] = EKF(f,h,F,H,Q,R,t,u,y,x0,P0)
%
% Author: Tamas Kis
% Last Update: 2021-08-18
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (function_handle) f(x,u,t) --> (n×1) discrete dynamics
%   h       - (function_handle) h(x,t)   --> (n×1) discrete measurement
%   F       - (function_handle) F(x,u,t) --> (n×n) dynamics Jacobian
%   H       - (function_handle) H(x,t)   --> (m×n) sensitivity matrix
%   Q       - (n×n double) process noise covariance (assumed constant)
%   R       - (m×m double) measurement noise covariance (assumed constant)
%   t       - (T×1 or 1×T double) time vector
%   u       - (p×T double) control input time history
%   y       - (m×T double) measurement time history
%   x0      - (n×1 double) initial state estimate
%   P0      - (n×n double) initial error covariance
%
% -------
% OUTPUT:
% -------
%   x       - (n×T double) a posteriori state estimates
%   P       - (n×n×T double) a posteriori error covariances
%   tsol    - (1×1 double) average time for one filter iteration
%   rank_Ob - (T×1 double) rank of the observability matrix
%   z_pre   - (m×T double) pre-fit measurement residuals
%   z_post  - (m×T double) post-fit measurement residuals
%
%==========================================================================
function [x,P,tsol,rank_Ob,z_pre,z_post] = EKF(f,h,F,H,Q,R,t,u,y,x0,P0)
    
    % number of sample times, state dimension, and measurement dimension
    T = length(t);
    n = length(x0);
    m = length(y(:,1));
    
    % preallocates arrays
    x = zeros(n,T);
    P = zeros(n,n,T);
    rank_Ob = zeros(T,1);
    z_pre = zeros(m,T);
    z_post = zeros(m,T);
    
    % assigns initial conditions
    x(:,1) = x0;
    P(:,:,1) = P0;
    
    % filtering
    tic;
    for k = 2:T      
        
        % evaluates dynamics Jacobian at previous state estimate
        F_prev = F(x(:,k-1),u(:,k-1),t(k-1));
        
        % TIME UPDATE (PREDICT STEP)
        x_pred = f(x(:,k-1),u(:,k-1),t(k-1));
        P_pred = F_prev*P(:,:,k-1)*F_prev'+Q;
        
        % evaluates sensitivity matrix at predicted state
        Hk = H(x_pred,t(k));
        
        % pre-fit measurement residual
        z_pre(:,k) = y(:,k)-h(x_pred,t(k));
        
        % Kalman gain
        K = P_pred*Hk'/(Hk*P_pred*Hk'+R);
        
        % MEASUREMENT UPDATE (UPDATE STEP)
        x(:,k) = x_pred+K*z_pre(:,k);
        P(:,:,k) = (eye(n)-K*Hk)*P_pred;
        
        % post-fit measurement residual
        z_post(:,k) = y(:,k)-h(x(:,k),t(k));
        
        % rank of the observability matrix
        rank_Ob(k) = rank(obsv(F(x(:,k),u(:,k),t(k)),H(x(:,k),t(k))));
        
    end
    tsol = toc/(T-1);
    
end