%==========================================================================
%
% EKF_sim  Simulation of an extended Kalman filter with pre-computed
% measurements and control inputs.
%
%   [x,P,tsol,z_pre,z_post,Fk,Hk] = EKF_sim(fd,hd,F,H,Q,R,[],y,x0,P0)
%   [x,P,tsol,z_pre,z_post,Fk,Hk] = EKF_sim(fd,hd,F,H,Q,R,u,y,x0,P0)
%   [__] = EKF_sim(__,wb)
%
% Author: Tamas Kis
% Last Update: 2022-03-31
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   fd      - (1×1 function_handle) discrete dynamics equation,
%             xₖ₊₁ = fd(xₖ,uₖ,k) (fd : ℝⁿ×ℝᵐ×ℤ → ℝⁿ)
%   hd      - (1×1 function_handle) discrete measurement equation,
%             yₖ = hd(xₖ,k) (hd : ℝⁿ×ℤ → ℝᵖ)
%   F       - (1×1 function_handle) discrete dynamics Jacobian, 
%             Fₖ = F(xₖ,uₖ,k) (F : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   H       - (1×1 function_handle) discrete measurement Jacobian,
%             Hₖ = H(xₖ,k) (H : ℝⁿ×ℤ → ℝᵖˣⁿ)
%   Q       - (1×1 function_handle) process noise covariance, 
%             Qₖ = Q(xₖ,uₖ,k) (Q : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   R       - (1×1 function_handle) measurement noise covariance, 
%             Rₖ = R(xₖ,k) (R : ℝⁿ×ℤ → ℝᵖˣᵖ)
%   u       - (m×(N-1) double) (OPTIONAL) control input history
%   y       - (p×N double) measurement history
%   x0      - (n×1 double) initial state estimate
%   P0      - (n×n double) initial error covariance
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
%==========================================================================
function [x,P,tsol,z_pre,z_post,Fk,Hk] = EKF_sim(fd,hd,F,H,Q,R,u,y,x0,...
    P0,wb)
    
    % initializes the waitbar
    if (nargin == 11)
        [wb,prop,display_waitbar] = initialize_waitbar(wb,...
            'Running extended Kalman filter...');
    else
        display_waitbar = false;
    end
    
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
            Hk(:,:,kk)] = EKF(x(:,kk-1),P(:,:,kk-1),u(:,kk-1),y(:,kk),...
            fd,hd,F,H,Q,R,k);

        % updates waitbar
        if display_waitbar, prop = update_waitbar(k,N,wb,prop); end
        
    end
    tsol = toc/N;

    % closes waitbar
    if display_waitbar, close(wb); end
    
end