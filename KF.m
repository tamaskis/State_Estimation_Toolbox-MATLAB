% KF  Kalman filter.
%
%   [mu,Sigma,Ob] = KF(f,g,A,C,Q,R,u,y,mu0,Sigma0)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-22
%
%=========================================================================%
%
% INPUTS:
%   A       function handle for system (plant) matrix (n x n)
%   B       function handle for input (control) matrix (n x m)
%   C       function handle for measurement matrix (k x n)
%   Q       process noise covariance (n x n)
%   R       measurement noise covariance (k x k)
%   u       control input (m x T matrix of m x 1 vectors)
%   y       measurement (k x (T-1) matrix of k x 1 vectors)
%   mu0     initial mean (n x 1)
%   Sigma0 	initial covariance (n x n)
%
% OUTPUTS:
%   mu      posterior mean (n x T matrix of n x 1 vectors)
%   Sigma   posterior covariance (n x n x T array of n x n matrices)
%   Ob      observability matrices (nk x n x T array of nk x n matrices)
%
%=========================================================================%
function [mu,Sigma] = KF(A,B,C,Q,R,u,y,mu0,Sigma0)
    
    % determines length of time vector
    T = length(u);
    
    % determines n
    n = length(mu0);
    
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
        At = A(mu(:,t-1),u(:,t-1),t);
        Bt = B(mu(:,t-1),u(:,t-1),t);
        mu_predict = At*mu(:,t-1)+Bt*u(:,t-1);
        Sigma_predict = At*Sigma(:,:,t-1)*At'+Q;
        
        % Kalman gain
        Ct = C(mu_predict,t);
        K = Sigma_predict*Ct'*inv(Ct*Sigma_predict*Ct'+R);
        
        % update step
        mu(:,t) = mu_predict+K*(y(:,t)-Ct*mu_predict);
        Sigma(:,:,t) = (I-K*Ct)*Sigma_predict;
        
        % observability matrix
        Ob(:,:,t) = obsv(A(mu(:,t),u(:,t),t),C(mu(:,t),t));
        
    end
    
end