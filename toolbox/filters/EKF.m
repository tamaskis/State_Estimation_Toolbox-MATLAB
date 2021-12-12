%==========================================================================
%
% EKF  Extended Kalman filter.
%
%   [x,P,tsol,rank_Ob,z_pre,z_post] = EKF(f,h,F,H,Q,R,t,u,y,x0,P0)
%   [x,P,tsol,rank_Ob,z_pre,z_post] = EKF(f,h,F,H,Q,R,t,u,y,x0,P0,wb)
%
% Author: Tamas Kis
% Last Update: 2021-12-12
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (1×1 function_handle) f(x,u,k) --> discrete-time nonlinear
%             dynamics (f:Rn×Rm×R->Rn)
%   h       - (1×1 function_handle) h(x,k) --> discrete-time nonlinear 
%             measurement (h:Rn×R->Rp)
%   F       - (1×1 function_handle) F(x,u,k) --> dynamics Jacobian 
%             (F:Rn×Rm×R->Rn×Rn)
%   H       - (1×1 function_handle) H(x,k) -->  measurement Jacobian 
%             (H:Rn×R->Rp×Rp)
%   Q       - (n×n double) process noise covariance (assumed constant)
%   R       - (p×p double) measurement noise covariance (assumed constant)
%   t       - (T×1 double) time vector
%   u       - (m×T double) control input time history
%   y       - (p×T double) measurement time history
%   x0      - (n×1 double) initial state estimate
%   P0      - (n×n double) initial error covariance
%   wb      - (OPTIONAL) (char or 1×1 logical) waitbar parameters
%               --> input as "true" if you want waitbar with default 
%                   message displayed
%               --> input as a char array storing a message if you want a
%                   custom message displayed on the waitbar
%
% -------
% OUTPUT:
% -------
%   x       - (n×T double) a posteriori state estimates
%   P       - (n×n×T double) a posteriori error covariances
%   tsol    - (1×1 double) average time for one filter iteration
%   rank_Ob - (T×1 double) rank of the observability matrix
%   z_pre   - (p×T double) pre-fit measurement residuals
%   z_post  - (p×T double) post-fit measurement residuals
%
%==========================================================================
function [x,P,tsol,rank_Ob,z_pre,z_post] = EKF(f,h,F,H,Q,R,t,u,y,x0,P0,wb)
    
    % -------------------
    % Setting up waitbar.
    % -------------------
    
    % initializes the waitbar if "wb" is input
    if nargin == 12
        display_waitbar = true;
        [wb,prop] = init_waitbar(wb);

    % defaults waitbar to off if "wb" not input
    else
        display_waitbar = false;
    end

    % -----------------------
    % Extended Kalman filter.
    % -----------------------

    % number of sample times, state dimension, and measurement dimension
    T = length(t);
    n = length(x0);
    p = length(y(:,1));
    
    % preallocates arrays
    x = zeros(n,T);
    P = zeros(n,n,T);
    rank_Ob = zeros(T,1);
    z_pre = zeros(p,T);
    z_post = zeros(p,T);
    
    % assigns initial conditions
    x(:,1) = x0;
    P(:,:,1) = P0;
    
    % filtering
    tic;
    for k = 2:T      
        
        % evaluates dynamics Jacobian at previous state estimate
        F_prev = F(x(:,k-1),u(:,k-1),k-2);
        
        % TIME UPDATE (PREDICT STEP)
        x_pred = f(x(:,k-1),u(:,k-1),k-2);
        P_pred = F_prev*P(:,:,k-1)*F_prev'+Q;
        
        % evaluates sensitivity matrix at predicted state
        Hk = H(x_pred,k);
        
        % pre-fit measurement residual
        z_pre(:,k) = y(:,k)-h(x_pred,k);
        
        % Kalman gain
        K = P_pred*Hk'/(Hk*P_pred*Hk'+R);
        
        % MEASUREMENT UPDATE (UPDATE STEP)
        x(:,k) = x_pred+K*z_pre(:,k);
        P(:,:,k) = (eye(n)-K*Hk)*P_pred;
        
        % post-fit measurement residual
        z_post(:,k) = y(:,k)-h(x(:,k),k);
        
        % rank of the observability matrix
        rank_Ob(k) = rank(obsv(F(x(:,k),u(:,k),k),H(x(:,k),k)));

        % updates waitbar
        if display_waitbar, prop = update_waitbar(k,T,wb,prop); end
        
    end
    tsol = toc/(T-1);

    % closes waitbar
    if display_waitbar, close(wb); end

    % -------------
    % Subfunctions.
    % -------------

    %----------------------------------------------------------------------
    % update_waitbar  Updates the waitbar.
    %----------------------------------------------------------------------
    %
    % INPUT:
    %   wb      - (1×1 Figure) waitbar
    %   i       - (1×1 double) current iteration
    %  	N       - (1×1 double) total number of iterations
    %   prop    - (1×1 double) cutoff proportion to trigger waitbar update
    %
    % OUTPUT:
    %   prop    - (1×1 double) cutoff proportion to trigger waitbar update
    %
    % NOTE:
    %   --> "prop" is an integer multiple of 0.1 so that the waitbar is
    %       only updated after every additional 10% of progress.
    %
    %----------------------------------------------------------------------
    function prop = update_waitbar(i,N,wb,prop)
        
        % only updates waitbar if current proportion exceeds cutoff prop.
        if i/N > prop
            
            % updates waitbar
            waitbar(i/N,wb);
            
            % updates cutoff proportion needed to trigger waitbar update
            prop = prop+0.1;
            
        end
        
    end
    
    %----------------------------------------------------------------------
    % init_waitbar  Initialize the waitbar.
    %----------------------------------------------------------------------
    %
    % INPUT:
    %   wb      - (OPTIONAL) (char or 1×1 logical) waitbar parameters
    %               --> input as "true" if you want waitbar with default 
    %                   message displayed
    %               --> input as a char array storing a message if you want
    %                   a custom message displayed on the waitbar
    %
    % OUTPUT:
    %   wb      - (1×1 Figure) waitbar
    %   prop    - (1×1 double) cutoff proportion to trigger waitbar update,
    %             initialized to 0.1
    %
    %----------------------------------------------------------------------
    function [wb,prop] = init_waitbar(wb)

        % initializes cutoff proporation to 0
        prop = 0;

        % initializes waitbar with default message
        if islogical(wb)
            wb = waitbar(0,'Running extended Kalman filter...');

        % initializes waitbar with custom messagae
        else
            wb = waitbar(0,wb);
        end

    end
    
end