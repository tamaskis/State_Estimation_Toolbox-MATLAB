%==========================================================================
%
% UKF_predict  UKF predict step (time update).
%
%   [x_pred,P_pred] = UKF_predict(x_prev,P_prev,u_prev,k,fd,Q)
%
% Author: Tamas Kis
% Last Update: 2022-03-28
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x_prev  - (n×1 double) state estimate at previous sample time
%   P_prev  - (n×n double) error covariance at previous sample time
%   u_prev  - (m×1 double) control input at previous sample time
%   k       - (1×1 double) current sample number
%   fd      - (1×1 function_handle) discrete nonlinear dynamics equation,
%             xₖ₊₁ = fd(xₖ,uₖ,k) (fd : ℝⁿ×ℝᵐ×ℤ → ℝⁿ)
%   Q       - (1×1 function_handle) Qₖ = Q(xₖ,uₖ,k) --> process noise 
%             covariance (Q : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%
% -------
% OUTPUT:
% -------
%   x_pred  - (n×1 double) a priori state estimate at current sample time
%   P_pred  - (n×n double) a priori error covariance at current sample time
%
%==========================================================================
function [x_pred,P_pred] = UKF_predict(x_prev,P_prev,u_prev,k,fd,Q)
    
    % state dimension
    n = length(x_prev);

    % sigma points from state estimate statistics at previous sample time
    [Chi,w] = UT(x_prev,P_prev);
    
    % passing sigma points through nonlinear dynamics
    for i = 1:(2*n+1)
        Chi(:,i) = fd(Chi(:,i),u_prev,k-1);
    end
    
    % a priori state estimate at current sample time
    [x_pred,P_tilde] = iUT(Chi,w);
    
    % a priori error covariance at current sample time
    P_pred = P_tilde+Q(x_prev,u_prev,k-1);
    
end