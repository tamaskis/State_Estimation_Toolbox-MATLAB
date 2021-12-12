%==========================================================================
%
% measurement_jacobian  Numerically approximates (to within double 
% precision using the complex-step approximation) the measurement Jacobian.
%
%   Hk = measurement_jacobian(h,xk,k)
%
% See also dynamics_jacobian.
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
%   h       - (function_handle) nonlinear measurement equation h(xk,uk,k)
%             (h:Rn,Rm,R->Rp)
%   xk      - (n×1 double) state estimate at kth sample time
%   k       - (1×1 double) sample #
%
% -------
% OUTPUT:
% -------
%   Hk      - (n×n double) measurement Jacobian at kth sample time
%
%==========================================================================
function Hk = measurement_jacobian(h,xk,k)
    Hk = ijacobian(@(xk)h(xk,k),xk);
end