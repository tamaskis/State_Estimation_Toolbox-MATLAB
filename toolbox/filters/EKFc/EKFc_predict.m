%==========================================================================
%
% EKFc_predict  EKF predict step (time update) for continuous-time systems.
%
%   [x_pred,P_pred,F_prev] = EKFc_predict(x_prev,P_prev,u_prev,f,A,Q,k,...
%       dt,t0,method)
%
% Author: Tamas Kis
% Last Update: 2022-03-31
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x_prev  - (n×1 double) state estimate at previous sample time
%   P_prev  - (n×n double) error covariance at previous sample time
%   u_prev  - (m×1 double) control input at previous sample time
%   f       - (1×1 function_handle) continuous dynamics equation,
%             dx/dt = f(x,u,t) (f : ℝⁿ×ℝᵐ×ℝ → ℝⁿ)
%   A       - (1×1 function_handle) continuous dynamics Jacobian, 
%             A(t) = A(x,u,t) (A : ℝⁿ×ℝᵐ×ℝ → ℝⁿˣⁿ)
%   Q       - (1×1 function_handle) process noise covariance, 
%             Qₖ = Q(xₖ,uₖ,k) (Q : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   k       - (1×1 double) current sample number
%   dt      - (1×1 double) time step, Δt
%   t0      - (1×1 double) initial time, t₀
%   method  - (char) integration method --> 'Euler', 'RK2', 'RK2 Heun', 
%             'RK2 Ralston', 'RK3', 'RK3 Heun', 'RK3 Ralston', 'SSPRK3',
%             'RK4', 'RK4 Ralston', 'RK4 3/8'
%
% -------
% OUTPUT:
% -------
%   x_pred  - (n×1 double) a priori state estimate at current sample time
%   P_pred  - (n×n double) a priori error covariance at current sample time
%   F_prev  - (n×n double) discrete dynamics Jacobian at previous sample 
%             time
%
% -----
% NOTE:
% -----
%   --> This function makes use of the fact that Fₖ₋₁ = Φ(tₖ,tₖ₋₁).
%
%==========================================================================
function [x_pred,P_pred,F_prev] = EKFc_predict(x_prev,P_prev,u_prev,f,A,...
    Q,k,dt,t0,method)
    
    % current sample time
    tk = k2t_num(k,dt,t0);
    
    % previous sample time
    tk_1 = tk-dt;
    
    % discrete dynamics Jacobian at previous sample time and a priori state 
    % estimate at current sample time
    [F_prev,x_pred] = Af2stm_num(A,f,x_prev,u_prev,tk_1,dt,method);
    
    % process noise covariance at previous sample time
    Q_prev = Q(x_prev,u_prev,k-1);
    
    % a priori error covariance at current sample time
    P_pred = F_prev*P_prev*F_prev.'+Q_prev;
    
end