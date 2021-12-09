%==========================================================================
%
% discretize_dynamics  Discretization of continuous-time nonlinear
% measurement.
%
%   fd = discretize_measurement(f,x,u,k)
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
function hd = discretize_measurement(h,t0,dt)
    hd = @(x,k) h(x,t0+k*dt);
end