% UKF  Unscented Kalman filter.
%
%   [mu,Sigma] = UKF(f,g,Q,R,u,y,mu0,Sigma0)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-22
%
%=========================================================================%
%
% INPUTS:
%   f       function handle for nonlinear dynamics equation (n x 1)
%   g       function handle for nonlinear measurement equation (k x 1)
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
%
%=========================================================================%
function [mu,Sigma] = UKF(f,g,Q,R,u,y,mu0,Sigma0)
    
    % determines length of time vector
    T = length(u);
    
    % determines n and k
    n = length(mu0);
    k = length(y(:,1));
    
    % preallocates arrays
    mu = zeros(n,T);
    Sigma = zeros(n,n,T);
    
    % assigns initial conditions
    mu(:,1) = mu0;
    Sigma(:,:,1) = Sigma0;
    
    % filtering
    for t = 2:T
        
        % preallocates vectors and matrices for update step
        y_bar = zeros(k,2*n+1);
        y_hat = zeros(k,1);
        Sigma_y = zeros(k,k);
        Sigma_xy = zeros(n,k);
        
        % prediction step
        [chi,w] = unscented_transform(mu(:,t-1),Sigma(:,:,t-1));
        for i = 1:(2*n+1)
            chi(:,i) = f(chi(:,i),u(:,t-1));
        end
        [mu_predict,Sigma_bar] = inverse_unscented_transform(chi,w);
        Sigma_predict = Sigma_bar+Q;
        
        % update step
        [chi,w] = unscented_transform(mu_predict,Sigma_predict);
        for i = 1:(2*n+1)
            y_bar(:,i) = g(chi(:,i));
            y_hat = y_hat+w(i)*y_bar(:,i);
        end
        for i = 1:(2*n+1)
            Sigma_y = Sigma_y+w(i)*(y_bar(:,i)-y_hat)*(y_bar(:,i)-...
                y_hat)';
            Sigma_xy = Sigma_xy+w(i)*(chi(:,i)-mu_predict)*(y_bar(:,i)-...
                y_hat)';
        end
        Sigma_y = Sigma_y+R;
        mu(:,t) = mu_predict+Sigma_xy*inv(Sigma_y)*(y(:,t)-y_hat);
        Sigma(:,:,t) = Sigma_predict-Sigma_xy*inv(Sigma_y)*Sigma_xy';

    end
    
end