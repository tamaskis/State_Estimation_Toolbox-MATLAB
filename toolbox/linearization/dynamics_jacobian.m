%==========================================================================
%
% dynamics_jacobian  Numerically approximates (to within double precision 
% using the complex-step approximation) the dynamics Jacobian.
%
%   Fk = dynamics_jacobian(fd,xk,uk,k)
%
% See also measurement_jacobian.
%
% Author: Tamas Kis
% Last Update: 2021-11-13
%
% REFERENCES:
%   [1] TODO
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   fd      - (function_handle) discrete-time nonlinear dynamics equation 
%             fd(xk,uk,k) (fd:Rn,Rm,R->Rn)
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
function Fk = dynamics_jacobian(fd,xk,uk,k)
    Fk = ijacobian(@(xk)fd(xk,uk,k),xk);
end