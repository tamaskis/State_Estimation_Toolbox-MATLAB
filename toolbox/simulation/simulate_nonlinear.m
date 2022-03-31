%==========================================================================
%
% simulate_nonlinear  Simulation of a nonlinear system.
%
%   [x,y] = simulate_nonlinear(fd,hd,Q,R,u,IC)
%   [x,y] = simulate_nonlinear(fd,hd,Q,R,u,IC,seed)
%
% Author: Tamas Kis
% Last Update: 2022-03-30
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
%   Q       - (1×1 function_handle) Qₖ = Q(xₖ,uₖ,k) --> process noise 
%             covariance (Q : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   R       - (1×1 function_handle) Rₖ = R(xₖ,k) --> measurement noise 
%             covariance (R : ℝⁿ×ℤ → ℝᵖˣᵖ)
%   u       - (m×(N-1) double) (OPTIONAL) control input history
%   IC      - (1×1 struct) structure storing initial conditions
%       • x0        - (n×1 double) initial state estimate
%       • P0        - (n×n double) initial error covariance
%                           OR
%       • x0_true   - (n×1 double) ground truth initial state
%   seed    - (1×1 double) (OPTIONAL) seed for random number generator
%
% -------
% OUTPUT:
% -------
%   x       - (n×N double) ground truth state trajectory
%   y       - (p×N double) ground truth measurement history
%
%==========================================================================
function [x,y] = simulate_nonlinear(fd,hd,Q,R,u,IC,seed)
    
    % sets random seed if specified
    if nargin == 8
        rng(seed);
    end
    
    % determines initial condition (if ground truth initial state not
    % given, the initial condition is sample from the initial state
    % statistics)
    if isfield(IC,'x0_true')
        x0_true = IC.x0_true;
    else
        x0_true = mvnrnd(IC.x0,IC.P0);
    end
    
    % number of samples, state dimension, and measurement dimension
    N = length(u);
    n = length(x0_true);
    p = size(R(x0_true,1),1);
    
    % preallocates arrays
    x = zeros(n,N);
    y = zeros(p,N);
    
    % assigns initial condition
    x(:,1) = x0_true;
    
    % simulation of ground truth
    for k = 2:N
        
        % state propagation with process noise
        x(:,k) = fd(x(:,k-1),u(:,k-1),k-1)+mvnrnd(zeros(n,1),Q)';
        
        % measurement simulation with measurement noise
        y(:,k) = hd(x(:,k),k)+mvnrnd(zeros(p,1),R)';
        
    end
    
end