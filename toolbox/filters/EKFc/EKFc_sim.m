%==========================================================================
%
% EKF_sim  Simulation of an extended Kalman filter with pre-computed
% measurements and control inputs.
%
%   [x,P,tsol,rank_Ob,z_pre,z_post] = EKF_sim(fd,hd,F,H,Q,R,[],y,x0,P0,dt)
%   [x,P,tsol,rank_Ob,z_pre,z_post] = EKF_sim(fd,hd,F,H,Q,R,u,y,x0,P0,dt)
%   [__] = EKF_sim(__,t0,method,wb)
%
% Author: Tamas Kis
% Last Update: 2022-03-31
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (1×1 function_handle) continuous dynamics equation,
%             dx/dt = f(x,u,t) (f : ℝⁿ×ℝᵐ×ℝ → ℝⁿ)
%   h       - (1×1 function_handle) continuous measurement equation,
%             y = h(x,t) (h : ℝⁿ×ℝ → ℝᵖ)
%   A       - (1×1 function_handle) continuous dynamics Jacobian, 
%             A(t) = A(x,u,t) (A : ℝⁿ×ℝᵐ×ℝ → ℝⁿˣⁿ)
%   C       - (1×1 function_handle) continuous measurement Jacobian,
%             C(t) = C(x,t) (C : ℝⁿ×ℝ → ℝᵖˣⁿ)
%   Q       - (1×1 function_handle) process noise covariance, 
%             Qₖ = Q(xₖ,uₖ,k) (Q : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   R       - (1×1 function_handle) measurement noise covariance, 
%             Rₖ = R(xₖ,k) (R : ℝⁿ×ℤ → ℝᵖˣᵖ)
%   u       - (m×(N-1) double) (OPTIONAL) control input history
%   y       - (p×N double) measurement history
%   x0      - (n×1 double) initial state estimate
%   P0      - (n×n double) initial error covariance
%   dt      - (1×1 double) time step, Δt
%   t0      - (1×1 double) (OPTIONAL) initial time, t₀
%   method  - (char) (OPTIONAL) integration method --> 'Euler', 'RK2', 
%             'RK2 Heun', 'RK2 Ralston', 'RK3', 'RK3 Heun', 'RK3 Ralston', 
%             'SSPRK3', 'RK4', 'RK4 Ralston', 'RK4 3/8' (defaults to 
%             'Euler')
%   wb      - (char or 1×1 logical) (OPTIONAL) waitbar parameters
%               --> input as "true" if you want waitbar with default 
%                   message displayed
%               --> input as a char array storing a message if you want a
%                   custom message displayed on the waitbar
%
% -------
% OUTPUT:
% -------
%   x       - (n×N double) a posteriori state estimates
%   P       - (n×n×N double) a posteriori error covariances
%   tsol    - (1×1 double) average time for one filter iteration
%   rank_Ob - (N×1 double) rank of the observability matrix
%   z_pre   - (p×N double) pre-fit measurement residuals
%   z_post  - (p×N double) post-fit measurement residuals
%
%==========================================================================
function [x,P,tsol,rank_Ob,z_pre,z_post] = EKFc_sim(f,h,A,C,Q,R,u,y,x0,...
    P0,dt,t0,method,wb)
    
    % -------------------------------
    % Defaulting optional parameters.
    % -------------------------------

    % defaults initial time to empty vector if not input
    if (nargin < 12)
        t0 = [];
    end

    % defaults integration method to empty vector if not input
    if (nargin < 13)
        method = [];
    end
    
    % -------------------
    % Setting up waitbar.
    % -------------------
    
    % determines if waitbar is on or off
    if (nargin < 14) || (islogical(wb) && ~wb)
        display_waitbar = false;
    else
        display_waitbar = true;
    end

    % sets the waitbar message (defaults to 'Running extended Kalman 
    % filter...')
    if display_waitbar
        if ischar(wb)
            msg = wb;
        else
            msg = 'Running extended Kalman filter...';
        end
    end

    % initialize cutoff proportion needed to trigger waitbar update to 0.1
    if display_waitbar, prop = 0.1; end

    % initializes the waitbar
    if display_waitbar, wb = waitbar(0,msg); end

    % -----------------------
    % Extended Kalman filter.
    % -----------------------
    
    % number of sample times (including initial time)
    N = size(y,2);

    % state (n) and measurement (p) dimensions
    n = length(x0);
    p = size(y,1);
    
    % defaults "u" to an array of zeros if not input
    if isempty(u)
        u = zeros(1,N-1);
    end
    
    % preallocates arrays
    x = zeros(n,N);
    P = zeros(n,n,N);
    rank_Ob = zeros(N,1);
    z_pre = zeros(p,N);
    z_post = zeros(p,N);
    Fk = zeros(n,n,N);
    Hk = zeros(p,n,N);
    
    % assigns initial conditions
    x(:,1) = x0;
    P(:,:,1) = P0;
    
    % filtering
    tic;
    for k = 1:(N-1)

        % index for accessing arrays (switch from 0- to 1-based indexing)
        kk = k+1;
        
        % one iteration of the extended Kalman filter
        [x(:,kk),P(:,:,kk),z_pre(:,kk),z_post(:,kk),Fk(:,:,kk-1),...
            Hk(:,:,kk)] = EKFc(x(:,kk-1),P(:,:,kk-1),u(:,kk-1),y(:,kk),...
            f,h,A,C,Q,R,k,dt,t0,method);

        % updates waitbar
        if display_waitbar, prop = update_waitbar(k,N,wb,prop); end
        
    end
    tsol = toc/N;

    % rank of the observability matrix
    for k = 0:(N-1)

        % index for accessing arrays (switch from 0- to 1-based indexing)
        kk = k+1;

        % observability at initial sample time
        if (k == 0)
            rank_Ob(kk) = rank(obsv(Fk(:,:,1),C(x0,0)));

        % observability at final sample time (set to NaN because control
        % input will not be known at final sample time)
        elseif (kk == N-1)
            rank_Ob(kk) = NaN;

        % observability at all other sample times
        else
            rank_Ob(kk) = rank(obsv(Fk(:,:,kk),Hk(:,:,kk)));

        end

    end

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
    %   n       - (1×1 double) current sample number (i.e. iteration)
    %  	N       - (1×1 double) total number of samples (i.e. iterations)
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
    
end