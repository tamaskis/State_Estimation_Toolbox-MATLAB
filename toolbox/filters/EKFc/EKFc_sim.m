%==========================================================================
%
% EKFc_sim  Simulation of an extended Kalman filter with pre-computed
% measurements and control inputs for continuous-time systems.
%
%   [x,P,tsol,z_pre,z_post] = EKFc_sim(fd,hd,F,H,Q,R,[],y,x0,P0,dt)
%   [x,P,tsol,z_pre,z_post] = EKFc_sim(fd,hd,F,H,Q,R,u,y,x0,P0,dt)
%   [__] = EKFc_sim(__,t0,method,wb)
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
%   z_pre   - (p×N double) pre-fit measurement residuals
%   z_post  - (p×N double) post-fit measurement residuals
%   Fk      - (n×n×N double) discrete dynamics Jacobians
%   Hk      - (p×p×N double) discrete measurement Jacobians
%
% -----
% NOTE:
% -----
%   --> N = number of samples
%   --> "Fk" and "Hk" are returned to aid observability analyses.
%
%==========================================================================
function [x,P,tsol,z_pre,z_post,Fk,Hk] = EKFc_sim(f,h,A,C,Q,R,u,y,x0,...
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
    
    % initializes the waitbar --> TODO: rewrite as script so all functions
    % within scope
    if (nargin == 14)
        [wb,prop,display_waitbar] = initialize_waitbar(wb,...
            'Running extended Kalman filter...');
    else
        display_waitbar = false;
    end

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

    % closes waitbar
    if display_waitbar, close(wb); end
    
end