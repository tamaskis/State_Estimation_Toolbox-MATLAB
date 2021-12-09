%==========================================================================
%
% KF  Kalman filter.
%
%   [x,P,tsol,rank_Ob,z_pre,z_post] = KF(Phi,Gamma,H,Q,R,t,u,y,x0,P0)
%   [x,P,tsol,rank_Ob,z_pre,z_post] = KF(Phi,Gamma,H,Q,R,t,u,y,x0,P0,wb)
%
% Author: Tamas Kis
% Last Update: 2021-11-15
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   Phi  	- (function_handle) Φ(x,u,k) --> (Φ:Rn×Rm×R->Rn×Rn) state 
%                                            transition matrix
%   Gamma   - (function_handle) Γ(x,u,k) --> (Γ:Rn×Rm×R->Rn×Rm) input
%                                            matrix
%   C       - (function_handle) C(x,k)   --> (C:Rn×R->Rp×Rn) measurement 
%                                            matrix
%   Q       - (n×n double) process noise covariance (assumed constant)
%   R       - (m×m double) measurement noise covariance (assumed constant)
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
function [x,P,tsol,rank_Ob,z_pre,z_post] = KF(Phi,Gamma,H,Q,R,t,u,y,x0,...
    P0,wb)
    
    % -------------------
    % Setting up waitbar.
    % -------------------
    
    % initializes the waitbar if "wb" is input
    if nargin == 11
        display_waitbar = true;
        [wb,prop] = init_waitbar(wb);

    % defaults waitbar to off if "wb" not input
    else
        display_waitbar = false;
    end

    % --------------
    % Kalman filter.
    % --------------

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
        
        % state transition and input matrices
        Phi_prev = Phi(x(:,k-1),u(:,k-1),k-1);
        Gamma_prev = Gamma(x(:,k-1),u(:,k-1),k-1);
        
        % TIME UPDATE (PREDICT STEP)
        x_pred = Phi_prev*x(:,k-1)+Gamma_prev*u(:,k-1);
        P_pred = Phi_prev*P(:,:,k-1)*Phi_prev'+Q;
        
        % measurement matrix
        Hk = H(x_pred,k);
        
        % pre-fit residual
        z_pre(:,k) = y(:,k)-Hk*x_pred;
        
        % Kalman gain
        K = P_pred*Hk'/(Hk*P_pred*Hk'+R);
        
        % MEASUREMENT UPDATE (UPDATE STEP)
        x(:,k) = x_pred+K*z_pre(:,k);
        P(:,:,k) = (eye(n)-K*Hk)*P_pred;
        
        % post-fit residual
        z_post(:,k) = y(:,k)-Hk*x(:,k);
        
        % rank of the observability matrix
        rank_Ob(k) = rank(obsv(Phi(x(:,k),u(:,k),k),H(x(:,k),k)));

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
            wb = waitbar(0,'Running Kalman filter...');

        % initializes waitbar with custom messagae
        else
            wb = waitbar(0,wb);
        end

    end
    
end