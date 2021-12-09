%==========================================================================
%
% discretize_dynamics  Discretization of continuous-time nonlinear
% dynamics.
%
%   fd = discretize_dynamics(f,x,u,k)
%
% Author: Tamas Kis
% Last Update: 2021-11-15
%
% REFERENCES:
%   [1] TODO
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (function_handle) continuous-time nonlinear dynamics equation 
%             f(x,u,t) (fd:Rn,Rm,R->Rn)
%   xk      - (n×1 double) state estimate at kth sample time
%   uk      - (m×1 double) control vector at kth sample time
%   k       - (1×1 double) sample #
%
% -------
% OUTPUT:
% -------
%   Fk      - (n×n double) dynamics Jacobian at kth sample time
%
%==========================================================================
function fd = discretize_dynamics(f,t0,dt)
    fd = @(x,u,k) x+f(x,u,t0+k*dt)*dt;
end