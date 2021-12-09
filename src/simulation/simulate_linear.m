%==========================================================================
%
% simulate_linear  Simulation of a linear system.
%
%   [x,y] = simulate_linear(Phi,Gamma,H,Q,R,u,IC)
%   [x,y] = simulate_linear(Phi,Gamma,H,Q,R,u,IC,seed)
%
% Author: Tamas Kis
% Last Update: 2021-11-14
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
%   R       - (p×p double) measurement noise covariance (assumed constant)
%   u       - (m×T double) control input time history
%   IC      - (struct) structure storing initial conditions
%       • x0        - (n×1 double) initial state estimate
%       • P0        - (n×n double) initial error covariance
%                           OR
%       • x0_true   - (n×1 double) ground truth initial state
%   seed    - (OPTIONAL) (1×1 double) seed for random number generator
%
% -------
% OUTPUT:
% -------
%   x       - (n×T double) ground truth dynamics simulation
%   y       - (p×T double) ground truth measurement simulation
%
%==========================================================================
function [x,y] = simulate_linear(Phi,Gamma,C,Q,R,u,IC,seed)
    
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
    p = size(R,1);
    
    % preallocates arrays
    x = zeros(n,T);
    y = zeros(p,T);
    
    % assigns initial condition
    x(:,1) = x0_true;
    
    % simulation of ground truth
    for k = 2:T
        
        % state transition, input, and sensitivity matrices
        Phi_prev = Phi(x(:,k-1),u(:,k-1),k-1);
        Gamma_prev = Gamma(x(:,k-1),u(:,k-1),k-1);
        Ck = C(x(:,k),k);
        
        % state propagation with process noise
        x(:,k) = Phi_prev*x(:,k-1)+Gamma_prev*u(:,k-1)+...
            mvnrnd(zeros(n,1),Q)';
        
        % measurement simulation with measurement noise
        y(:,k) = Ck*x(:,k)+mvnrnd(zeros(p,1),R)';
        
    end
    
end