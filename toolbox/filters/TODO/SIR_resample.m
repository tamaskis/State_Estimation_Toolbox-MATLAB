% PF  Sequential importance resampling (SIR) particle filter.
%
%   [mu,Sigma] = PF(f,g,w,pv,u,y,x0,N)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-22
%
%=========================================================================%
%
% INPUTS:
%   f       function handle for nonlinear dynamics equation (n x 1)
%   g       function handle for nonlinear measurement equation (k x 1)
%   w       function handle for randomly sampling process noise (n x 1)
%   pv      function handle for measurement noise PDF (k x 1)
%   u       control input (m x T matrix of m x 1 vectors)
%   y       measurement (k x (T-1) matrix of k x 1 vectors)
%   mu0     initial mean (n x 1)
%   Sigma0 	initial covariance (n x n)
%   N       number of particles
%
% OUTPUTS:
%   mu      posterior mean (n x T matrix of n x 1 vectors)
%   Sigma   posterior covariance (n x n x T array of n x n matrices)
%
%=========================================================================%
function [mu,Sigma] = SIR(f,g,w,pv,u,y,x0,N)
    
    % determines length of time vector
    T = length(u);
    
    % determines n
    n = size(x0,1);
    
    % preallocates arrays
    mu = zeros(n,T);
    Sigma = zeros(n,n,T);
    x_predict = zeros(n,N);
    w_bar = zeros(N,1);
    w_hat = zeros(N,1);
    
    % initial particle set
    x = x0;
    
    % filtering
    for t = 2:T
        
        % prediction step
        for i = 1:N
            x_predict(:,i) = f(x(:,i),u(:,t-1))+w();
        end
        
        % update step
        for i = 1:N
            w_hat(i) = pv(y(:,t)-g(x_predict(:,i)));
        end
        w_sum = sum(w_hat);
        for i = 1:N
            w_bar(i) = w_hat(i)/w_sum;
        end
        
        % resample step
        for j = 1:N
            s = rand;
            w_cumulative_sum = 0;
            i = 1;
            new_particle = false;
            while (i < N) && (new_particle == false)
                w_cumulative_sum = w_cumulative_sum+w_bar(i);
                if s < w_cumulative_sum
                    x(:,j) = x_predict(:,i);
                    new_particle = true;
                end
                i = i+1;
            end
        end
        
        % resets the weights to 1/N
        w_bar(:) = 1/N;
        
        % posterior (empirical) mean
        for i = 1:N
            mu(:,t) = mu(:,t)+w_bar(i)*x(:,i);
        end
        
        % posterior (empirical) covariance
        for i = 1:N
            Sigma(:,:,t) = Sigma(:,:,t)+w_bar(i)*(x(:,i)-mu(:,t))*...
                (x(:,i)-mu(:,t))';
        end
        
    end
    
end