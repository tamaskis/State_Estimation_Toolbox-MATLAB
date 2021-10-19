%==========================================================================
%
% UKF  Unscented Kalman filter.
%
%   [x,P,tsol,z_pre,z_post] = UKF(f,h,Q,R,t,u,y,x0,P0)
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
%   h       - (function_handle) h(x,t)   --> (m×1) discrete measurement
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
%   z_pre   - (m×T double) pre-fit measurement residuals
%   z_post  - (m×T double) post-fit measurement residuals
%
%==========================================================================
function [x,P,tsol,z_pre,z_post] = UKF(f,h,Q,R,t,u,y,x0,P0)
    
    % number of sample times, state dimension, and measurement dimension
    T = length(t);
    n = length(x0);
    m = length(y(:,1));
    
    % preallocates arrays
    x = zeros(n,T);
    P = zeros(n,n,T);
    z_pre = zeros(m,T);
    z_post = zeros(m,T);
    
    % assigns initial conditions
    x(:,1) = x0;
    P(:,:,1) = P0;
    
    % filtering
    tic;
    for k = 2:T
        
        % sigma points from previous state estimate statistics
        [Chi,w] = unscented_transform(x(:,k-1),P(:,:,k-1));
        
        % passing sigma points through nonlinear dynamics
        for i = 1:(2*n+1)
            Chi(:,i) = f(Chi(:,i),u(:,k-1),t(k-1));
        end
        
        % TIME UPDATE (PREDICT STEP)
        [x_pred,P_tilde] = inverse_unscented_transform(Chi,w);
        P_pred = P_tilde+Q;
        
        % sigma points from predicted state estimate statistics
        [Chi,w] = unscented_transform(x_pred,P_pred);
        
        % passing sigma points through nonlinear measurement equation
        Y = zeros(m,2*n+1);
        y_pred = zeros(m,1);
        for i = 1:(2*n+1)
            Y(:,i) = h(Chi(:,i),t(k));
            y_pred = y_pred+w(i)*Y(:,i);
        end
        
        % pre-fit residual
        z_pre(:,k) = y(:,k)-y_pred;
        
        % covariance of predicted measurement and cross covariance between
        % predicted state and predicted measurement
        Py = zeros(m,m);
        Pxy = zeros(n,m);
        for i = 1:(2*n+1)
            Py = Py+w(i)*(Y(:,i)-y_pred)*(Y(:,i)-y_pred)';
            Pxy = Pxy+w(i)*(Chi(:,i)-x_pred)*(Y(:,i)-y_pred)';
        end
        Py = Py+R;
        
        % Kalman gain
        K = Pxy/Py;
        
        % MEASUREMENT UPDATE (UPDATE STEP)
        x(:,k) = x_pred+K*z_pre(:,k);
        P(:,:,k) = P_pred-K*Py*K';
        
        % post-fit residual
        z_post(:,k) = y(:,k)-h(x(:,k),t(k));

    end
    tsol = toc/(T-1);
    
end