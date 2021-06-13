% EKS_angle  Extended Kalman smoother for angular measurements.
%
%   [mu,Sigma] = EKS_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,i_angle)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-22
%
%=========================================================================%
%
% INPUTS:
%   f           function handle for nonlinear dynamics equation (n x 1)
%   g           function handle for nonlinear measurement equation (k x 1)
%   A           function handle for dynamics Jacobian (n x n)
%   C           function handle for measurement Jacobian (k x n)
%   Q           process noise covariance (n x n)
%   R           measurement noise covariance (k x k)
%   u           control input (m x T matrix of m x 1 vectors)
%   y           measurement (k x (T-1) matrix of k x 1 vectors)
%   mu0         initial mean (n x 1)
%   Sigma0      initial covariance (n x n)
%   i_angle     vector of indices corresponding to the elements of the
%               measurement vector that store angular measurements
%
% OUTPUTS:
%   mu          posterior mean (n x T matrix of n x 1 vectors)
%   Sigma       posterior covariance (n x n x T array of n x n matrices)
%
%=========================================================================%
function [mu,Sigma] = EKS_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,i_angle)
    
    % forward pass
    [mu_forward,Sigma_forward] = EKF_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,...
        i_angle);
    
    % determines length of time vector
    T = length(u);
    
    % preallocates arrays for backward pass
    mu = zeros(size(mu_forward));
    Sigma = zeros(size(Sigma_forward));
    
    % sets "initial conditions" for mu and Sigma at time step T
    mu(:,T) = mu_forward(:,T);
    Sigma(:,:,T) = Sigma_forward(:,:,T);
    
    % performs backward pass
    for t = (T-1):(-1):1
        
        % dynamics Jacobian based on forward pass
        At = A(mu_forward(:,t),u(:,t),t);
        
        % predicted covariance at time t+1
        Sigma_t1 = At*Sigma_forward(:,:,t)*At'+Q;
        
        % "backward" Kalman gain
        K = Sigma_forward(:,:,t)*At'*inv(Sigma_t1);
        
        % corrected mean
        mu(:,t) = mu_forward(:,t)+K*(mu(:,t+1)-f(mu_forward(:,t),u(:,t),...
            t));
        
        % corrected covariance
        Sigma(:,:,t) = Sigma_forward(:,:,t)+K*(Sigma(:,:,t+1)-Sigma_t1)*K';
        
    end
    
end