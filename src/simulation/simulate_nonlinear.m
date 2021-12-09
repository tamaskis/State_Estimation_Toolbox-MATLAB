%==========================================================================
%
% simulate_nonlinear  Simulation of a nonlinear system.
%
%   [x,y] = simulate_nonlinear(f,h,Q,R,t,u,IC)
%   [x,y] = simulate_nonlinear(f,h,Q,R,t,u,IC,seed)
%
% Author: Tamas Kis
% Last Update: 2021-11-14
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (function_handle) f(x,u,k) --> (f:Rn×Rm×R->Rn) discrete-time 
%                                            nonlinear dynamics equation
%   h       - (function_handle) h(x,k)   --> (h:Rn×R->Rp) discrete-time 
%                                            nonlinear measurement equation
%   Q       - (n×n double) process noise covariance (assumed constant)
%   R       - (p×p double) measurement noise covariance (assumed constant)
%   t       - (T×1 double) time vector
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
function [x,y] = simulate_nonlinear(f,h,Q,R,t,u,IC,seed)
    
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
    
    % length of time vector, state dimension, and measurement dimension
    T = length(t);
    n = length(x0_true);
    p = size(R,1);
    
    % preallocates arrays
    x = zeros(n,T);
    y = zeros(p,T);
    
    % assigns initial condition
    x(:,1) = x0_true;
    
    % simulation of ground truth
    for k = 2:T
        
        % state propagation with process noise
        x(:,k) = f(x(:,k-1),u(:,k-1),k-1)+mvnrnd(zeros(n,1),Q)';
        
        % measurement simulation with measurement noise
        y(:,k) = h(x(:,k),k)+mvnrnd(zeros(p,1),R)';
        
    end
    
end