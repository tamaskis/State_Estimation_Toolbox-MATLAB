%==========================================================================
%
% EKFc  Extended Kalman filter (single iteration) for continuous-time 
% systems.
%
%   [xk,Pk,z_pre,z_post,F_prev,Hk] = EKFc(x_prev,P_prev,u_prev,yk,f,h,A,...
%       C,Q,R,k,dt,t0,method)
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
%   yk      - (p×1 double) measurement at current sample time
%   f       - (1×1 function_handle) continuous dynamics equation,
%             dx/dt = f(x,u,t) (f : ℝⁿ×ℝᵐ×ℝ → ℝⁿ)
%   h       - (1×1 function_handle) continuous measurement equation,
%             y = h(x,t) (h : ℝⁿ×ℝ → ℝᵖ)
%   A       - (1×1 function_handle) continuous dynamics Jacobian, 
%             A(t) = A(x,u,t) (A : ℝⁿ×ℝᵐ×ℝ → ℝⁿˣⁿ)
%   C       - (1×1 function_handle) continuous measurement Jacobian,
%             C(t) = C(x,t) (C : ℝⁿ×ℝ → ℝᵖˣⁿ)
%   Q       - (1×1 function_handle) process noise covariance, 
%             Qₖ = Q(xₖ,uₖ,k) (Q : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   R       - (1×1 function_handle) measurement noise covariance, 
%             Rₖ = R(xₖ,k) (R : ℝⁿ×ℤ → ℝᵖˣᵖ)
%   k       - (1×1 double) current sample number
%   dt      - (1×1 double) time step, Δt
%   t0      - (1×1 double) initial time, t₀
%   method  - (char) integration method --> 'Euler', 'RK2', 'RK2 Heun', 
%             'RK2 Ralston', 'RK3', 'RK3 Heun', 'RK3 Ralston', 'SSPRK3',
%             'RK4', 'RK4 Ralston', 'RK4 3/8' (defaults to 'Euler')
%
% -------
% OUTPUT:
% -------
%   xk      - (n×1 double) a posteriori state estimate
%   Pk      - (n×n double) a posteriori error covariance
%   z_pre   - (p×1 double) pre-fit measurement residual
%   z_post  - (p×1 double) post-fit measurement residual
%   F_prev  - (n×n double) discrete dynamics Jacobian at previous sample 
%             time
%   Hk      - (p×m double) discrete measurement Jacobian at current sample
%             time
%
%==========================================================================
function [xk,Pk,z_pre,z_post,F_prev,Hk] = EKFc(x_prev,P_prev,u_prev,yk,...
    f,h,A,C,Q,R,k,dt,t0,method)
    
    % predict step (time update)
    [x_pred,P_pred,F_prev] = EKFc_predict(x_prev,P_prev,u_prev,f,A,Q,k,...
        dt,t0,method);

    % update step (measurement update)
    [xk,Pk,z_pre,z_post,Hk] = EKFc_update(x_pred,P_pred,yk,h,C,R,k,dt,t0);
    
end