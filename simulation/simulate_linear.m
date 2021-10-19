%==========================================================================
%
% simulate_linear  Simulation of a linear system.
%
%   [x,y] = simulate_linear(Phi,Gamma,H,Q,R,u,IC)
%   [x,y] = simulate_linear(Phi,Gamma,H,Q,R,u,IC,seed)
%
% Last Update: 2021-07-27
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   Phi  	- (function_handle) Φ(x,u) --> (n×n) state transition matrix
%   Gamma   - (function_handle) Γ(x,u) --> (n×p) input matrix
%   H       - (function_handle) H(x)   --> (m×n) sensitivity matrix
%   Q       - (n×n) process noise covariance (assumed constant)
%   R       - (m×m) measurement noise covariance (asssumed constant)
%   u       - (p×T) control input time history
%   IC      - (struct) structure storing initial conditions
%       x0     	- (n×1) initial state estimate
%       P0      - (n×n) initial error covariance
%                           OR
%       x0_true - (n×1) ground truth initial state
%   seed    - (OPTIONAL) (1×1) seed for random number generator
%
% -------
% OUTPUT:
% -------
%   x       - (n×T) ground truth dynamics simulation
%   y       - (m×T) ground truth measurement simulation
%
%==========================================================================
function [x,y] = simulate_linear(Phi,Gamma,H,Q,R,u,IC,seed)
    
    % seeds random number generator
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
    
    % length of time vector, state dimension, and measurement dimension
    T = length(u);
    n = length(x0_true);
    m = size(R,1);
    
    % preallocates arrays
    x = zeros(n,T);
    y = zeros(m,T);
    
    % assigns initial condition
    x(:,1) = x0_true;
    
    % simulation of ground truth
    for k = 2:T
        
        % state transition, input, and sensitivity matrices
        Phi_prev = Phi(x(:,k-1),u(:,k-1));
        Gamma_prev = Gamma(x(:,k-1),u(:,k-1));
        Hk = H(x(:,k));
        
        % state propagation with process noise
        x(:,k) = Phi_prev*x(:,k-1)+Gamma_prev*u(:,k-1)+...
            mvnrnd(zeros(n,1),Q)';
        
        % measurement simulation with measurement noise
        y(:,k) = Hk*x(:,k)+mvnrnd(zeros(m,1),R)';
        
    end
    
end