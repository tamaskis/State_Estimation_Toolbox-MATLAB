%==========================================================================
%
% EKF  Extended Kalman filter.
%
%   [x,P,tsol,rank_Ob,z_pre,z_post] = EKF(f,h,F,H,Q,R,t,u,y,x0,P0)
%   [x,P,tsol,rank_Ob,z_pre,z_post] = EKF(f,h,F,H,Q,R,t,u,y,x0,P0,wb)
%
% Author: Tamas Kis
% Last Update: 2021-12-12
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (1×1 function_handle) f(x,u,k) --> discrete-time nonlinear
%             dynamics (f:Rn×Rm×R->Rn)
%   h       - (1×1 function_handle) h(x,k) --> discrete-time nonlinear 
%             measurement (h:Rn×R->Rp)
%   F       - (1×1 function_handle) F(x,u,k) --> dynamics Jacobian 
%             (F:Rn×Rm×R->Rn×Rn)
%   H       - (1×1 function_handle) H(x,k) -->  measurement Jacobian 
%             (H:Rn×R->Rp×Rp)
%   Q       - (n×n double) process noise covariance (assumed constant)
%   R       - (p×p double) measurement noise covariance (assumed constant)
%   t       - (T×1 double) time vector
%   u       - (m×T double) control input time history
%   y       - (p×T double) measurement time history
%   x0      - (n×1 double) initial state estimate
%   P0      - (n×n double) initial error covariance
%   wb      - (OPTIONAL) (char or 1×1 logical) waitbar parameters
%               --> input as "true" if you want waitbar with default 
%                   message displayed
%               --> input as a char array storing a message if you want a
%                   custom message displayed on the waitbar
%
% -------
% OUTPUT:
% -------
%   x       - (n×T double) a posteriori state estimates
%   P       - (n×n×T double) a posteriori error covariances
%   tsol    - (1×1 double) average time for one filter iteration
%   rank_Ob - (T×1 double) rank of the observability matrix
%   z_pre   - (p×T double) pre-fit measurement residuals
%   z_post  - (p×T double) post-fit measurement residuals
%
%==========================================================================
function [x,P,tsol,rank_Ob,z_pre,z_post] = EKF(f,h,F,H,Q,R,t,u,y,x0,P0,wb)
    
    

end