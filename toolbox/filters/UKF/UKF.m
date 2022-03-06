%==========================================================================
%
% UKF  Unscented Kalman filter.
%
%   [x,P,tsol,z_pre,z_post] = UKF(f,h,Q,R,t,u,y,x0,P0)
%   [x,P,tsol,z_pre,z_post] = UKF(f,h,Q,R,t,u,y,x0,P0,wb)
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
%             measurement equation (h:Rn×R->Rp)
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
%   z_pre   - (p×T double) pre-fit measurement residuals
%   z_post  - (p×T double) post-fit measurement residuals
%
%==========================================================================
function [x,P,tsol,z_pre,z_post] = UKF(f,h,Q,R,t,u,y,x0,P0,wb)
    
    % -------------------
    % Setting up waitbar.
    % -------------------
    
    % initializes the waitbar if "wb" is input
    if nargin == 10
        display_waitbar = true;
        [wb,prop] = init_waitbar(wb);

    % defaults waitbar to off if "wb" not input
    else
        display_waitbar = false;
    end

    % ------------------------
    % Unscented Kalman filter.
    % ------------------------

    % number of sample times, state dimension, and measurement dimension
    T = length(t);
    n = length(x0);
    p = length(y(:,1));
    
    % preallocates arrays
    x = zeros(n,T);
    P = zeros(n,n,T);
    z_pre = zeros(p,T);
    z_post = zeros(p,T);
    
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
            Chi(:,i) = f(Chi(:,i),u(:,k-1),k-1);
        end
        
        % TIME UPDATE (PREDICT STEP)
        [x_pred,P_tilde] = inverse_unscented_transform(Chi,w);
        P_pred = P_tilde+Q;
        
        % sigma points from predicted state estimate statistics
        [Chi,w] = unscented_transform(x_pred,P_pred);
        
        % passing sigma points through nonlinear measurement equation
        Y = zeros(p,2*n+1);
        y_pred = zeros(p,1);
        for i = 1:(2*n+1)
            Y(:,i) = h(Chi(:,i),k);
            y_pred = y_pred+w(i)*Y(:,i);
        end
        
        % pre-fit residual
        z_pre(:,k) = y(:,k)-y_pred;
        
        % covariance of predicted measurement and cross covariance between
        % predicted state and predicted measurement
        Py = zeros(p,p);
        Pxy = zeros(n,p);
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
        z_post(:,k) = y(:,k)-h(x(:,k),k);
        
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
            wb = waitbar(0,'Running unscented Kalman filter...');

        % initializes waitbar with custom messagae
        else
            wb = waitbar(0,wb);
        end

    end
    
end