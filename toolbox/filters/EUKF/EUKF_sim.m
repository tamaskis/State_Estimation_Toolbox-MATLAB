%==========================================================================
%
% EUKF_sim  Simulation of an extended/unscented Kalman filter with 
% pre-computed measurements and control inputs.
%
%   [x,P,tsol,z_pre,z_post] = EUKF_sim(fd,hd,F,Q,R,[],y,x0,P0)
%   [x,P,tsol,z_pre,z_post] = EUKF_sim(fd,hd,F,Q,R,u,y,x0,P0)
%   [__] = EUKF_sim(__,wb)
%
% Author: Tamas Kis
% Last Update: 2022-03-20
%
% REFERENCES:
%   [1] TODO
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   fd      - (1×1 function_handle) discrete nonlinear dynamics equation,
%             xₖ₊₁ = fd(xₖ,uₖ,k) (fd : ℝⁿ×ℝᵐ×ℤ → ℝⁿ)
%   hd      - (1×1 function_handle) discrete nonlinear measurement 
%             equation, yₖ = hd(xₖ,k) (fd : ℝⁿ×ℤ → ℝᵖ)
%   F       - (1×1 function_handle) Fₖ = F(xₖ,uₖ,k) --> discrete dynamics
%             Jacobian (F : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   Q       - (n×n double) process noise covariance (assumed constant)
%   R       - (p×p double) measurement noise covariance (assumed constant)
%   u       - (m×(N-1) double) (OPTIONAL) control input time history
%   y       - (p×N double) measurement time history
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
%
%==========================================================================
function [x,P,tsol,z_pre,z_post] = EUKF_sim(fd,hd,F,Q,R,u,y,x0,P0,wb)
    
    % -------------------
    % Setting up waitbar.
    % -------------------
    
    % determines if waitbar is on or off
    if (nargin < 10) || (islogical(wb) && ~wb)
        display_waitbar = false;
    else
        display_waitbar = true;
    end

    % sets the waitbar message (defaults to 'Running extended/unscented 
    % Kalman filter...')
    if display_waitbar
        if ischar(wb)
            msg = wb;
        else
            msg = 'Running extended/unscented Kalman filter...';
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
    z_pre = zeros(p,N);
    z_post = zeros(p,N);
    Fk = zeros(n,n,N);
    
    % assigns initial conditions
    x(:,1) = x0;
    P(:,:,1) = P0;
    
    % filtering
    tic;
    for k = 1:(N-1)

        % index for accessing arrays (switch from 0- to 1-based indexing)
        kk = k+1;
        
        % one iteration of the extended/unscented Kalman filter
        [x(:,kk),P(:,:,kk),z_pre(:,kk),z_post(:,kk)] = EUKF(x(:,kk-1),...
            P(:,:,kk-1),u(:,kk-1),y(:,kk),k,fd,hd,F,Q,R);

        % updates waitbar
        if display_waitbar, prop = update_waitbar(k,N,wb,prop); end
        
    end
    tsol = toc/N;
    
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