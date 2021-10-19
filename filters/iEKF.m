%==========================================================================
%
% iEKF  Iterated extended Kalman filter.
%
%   [x,P,tsol,rank_Ob,z_pre,z_post] = iEKF(f,h,F,H,Q,R,u,y,x0,P0)
%
% Author: Tamas Kis
% Last Update: 2021-08-18
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (function_handle) f(x,u) --> (n×1) discrete dynamics
%   h       - (function_handle) h(x)   --> (n×1) discrete measurement
%   F       - (function_handle) F(x,u) --> (n×n) dynamics Jacobian
%   H       - (function_handle) H(x)   --> (m×n) sensitivity matrix
%   Q       - (n×n) process noise covariance (assumed constant)
%   R       - (m×m) measurement noise covariance (assumed constant)
%   u       - (p×T) control input time history
%   y       - (m×T) measurement time history
%   x0      - (n×1) initial state estimate
%   P0      - (n×n) initial error covariance
%
% -------
% OUTPUT:
% -------
%   x       - (n×T) a posteriori state estimates
%   P       - (n×n×T) a posteriori error covariances
%   tsol    - (1×1) average time for one filter iteration
%   rank_Ob - (T×1) rank of the observability matrix
%   z_pre   - (m×T) pre-fit measurement residuals
%   z_post  - (m×T) post-fit measurement residuals
%
%==========================================================================
function [x,P,tsol,rank_Ob,z_pre,z_post] = iEKF(f,h,F,H,Q,R,u,y,x0,P0)
    
    % number of sample times, state dimension, and measurement dimension
    T = length(u);
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
        F_prev = F(x(:,k-1),u(:,k-1));
        
        % TIME UPDATE (PREDICT STEP)
        x_pred = f(x(:,k-1),u(:,k-1));
        P_pred = F_prev*P(:,:,k-1)*F_prev'+Q;
        
        % pre-fit measurement residual
        z_pre(:,k) = y(:,k)-h(x_pred);
        
        % initializes "old" and "new" estimates for iterated update step
        x_old = x_pred;
        x_new = 0;
        
        % sets error so loop is entered
        err = 1;
        
        % iteration
        i = 1;
        while (err > 1e-10) && (i < 1e3)
            
            % evaluates sensitivity matrix at predicted state
            Hi = H(x_old);
            
            % Kalman gain
            K = P_pred*Hi'/(Hi*P_pred*Hi'+R);
            
            % new predicted state
            x_new = x_pred+K*(y(:,k)-h(x_old))+K*Hi*(x_old-x_pred);
            
            % relative error
            err = norm(x_new-x_old);
            
            % stores new state estimate for next iteration
            x_old = x_new;
            
            % increments loop index
            i = i+1;
            
        end
        
        % MEASUREMENT UPDATE (UPDATE STEP)
        x(:,k) = x_new;
        P(:,:,k) = (eye(n)-K*Hi)*P_pred;
        
        % post-fit measurement residual
        z_post(:,k) = y(:,k)-h(x(:,k));
        
        % rank of the observability matrix
        rank_Ob(k) = rank(obsv(F(x(:,k),u(:,k)),H(x(:,k))));
        
    end
    tsol = toc/(T-1);
    
end