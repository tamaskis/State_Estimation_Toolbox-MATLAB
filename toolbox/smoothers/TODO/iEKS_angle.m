% iEKS_angle  Iterated extended Kalman smoother for angular measurements.
%
%   [mu,Sigma] = iEKS_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,i_angle)
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
function [mu,Sigma] = iEKS_angle(f,g,A,C,Q,R,u,y,mu0,Sigma0,i_angle)
    
    % initializes mean and covariance for iterative solution
    mu0_old = mu0;
    Sigma0_old = Sigma0;
        
    % sets error so loop is entered
    err = 1;
        
    % iterated extended Kalman smoother
    i = 1;
    while (err > 1e-10) && (i < 100)
        
        % runs iteration of EKS
        [mu,Sigma] = EKS_angle(f,g,A,C,Q,R,u,y,mu0_old,Sigma0_old,i_angle);
        
        % new smoothed prior distribution
        mu0_new = mu(:,1);
        Sigma0_new = Sigma(:,:,1);
        
        % error between iterations
        err = norm(mu0_new-mu0_old);
        
        % stores new prior distribution for next iteration
        mu0_old = mu0_new;
        Sigma0_old = Sigma0_new;
        
        % increments loop index
        i = i+1;
        
    end
    
end