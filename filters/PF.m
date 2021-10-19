%==========================================================================
%
% PF  Particle filter.
%
%   [mu,Sigma] = PF(f,g,w,pv,u,y,x0,N)
%
% Author: Tamas Kis
% Last Update: 2021-08-18
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (function_handle) f(x,u) --> (n×1) discrete dynamics
%   h       - (function_handle) h(x)   --> (n×1) discrete measurement
%   w       function handle for randomly sampling process noise (n x 1)
%   pv      function handle for measurement noise PDF (k x 1)
%   u       - (p×T) control input time history
%   y       - (m×T) measurement time history
%   x0      - (n×1) initial state estimate
%   P0      - (n×n) initial error covariance
%   N       - number of particles
%
% -------
% OUTPUT:
% -------
%   x       - (n×T) a posteriori state estimates
%   P       - (n×n×T) a posteriori error covariances
%   tsol    - (1×1) average time for one filter iteration
%   z_pre   - (m×T) pre-fit measurement residuals
%   z_post  - (m×T) post-fit measurement residuals
%
%==========================================================================
function [mu,Sigma] = PF(f,g,w,pv,u,y,x0,N)
    
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