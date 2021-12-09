% EKF_angle  Extended Kalman filter for angular measurements.
%
%   [mu,Sigma,Ob] = EKF_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,i_angle)
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
%   Ob          observability matrices (nk x n x T array of nk x n 
%               matrices)
%
%=========================================================================%
function [mu,Sigma,Ob] = EKF_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,i_angle)
    
    % determines length of time vector
    T = length(u);
    
    % determines n and k
    n = length(mu0);
    k = length(y(:,1));
    
    % identity matrix for update step
    I = eye(n);
    
    % preallocates arrays
    mu = zeros(n,T);
    Sigma = zeros(n,n,T);
    Ob = zeros(n*k,n,T);
    
    % assigns initial conditions
    mu(:,1) = mu0;
    Sigma(:,:,1) = Sigma0;
    
    % filtering
    for t = 2:T

        % prediction step
        mu_predict = f(mu(:,t-1),u(:,t-1),t-1);
        At = A(mu(:,t-1),u(:,t-1),t-1);
        Sigma_predict = At*Sigma(:,:,t-1)*At'+Q;
        
        % Kalman gain
        Ct = C(mu_predict,t);
        K = Sigma_predict*Ct'*inv(Ct*Sigma_predict*Ct'+R);
        
        % measurement residual
        z = y(:,t)-g(mu_predict,t);
        z(i_angle) = wrapToPi(z(i_angle));
        
        % update step
        mu(:,t) = mu_predict+K*z;
        Sigma(:,:,t) = (I-K*Ct)*Sigma_predict;
        
        % observability matrix
        Ob(:,:,t) = obsv(A(mu(:,t),u(:,t),t),C(mu(:,t),t));
        
    end
    
end