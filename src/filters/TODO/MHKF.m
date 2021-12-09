%==========================================================================
%
% MHKF  Multi-hypothesis Kalman filter.
%
%   [mu,Sigma] = MHKF(A,B,C,Q,R,u,y,mu0,Sigma0,Z)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
% INPUTS:
%   A       function handle for system (plant) matrix (n x n)
%   B       function handle for input (control) matrix (n x m)
%   C       function handle for measurement matrix (k x n)
%   Q       function handles for process noise covariance (n x n)
%   R       function haandle for measurement noise covariance (k x k)
%   u       control input (m x T matrix of m x 1 vectors)
%   y       measurement (k x (T-1) matrix of k x 1 vectors)
%   mu0     initial prior mean (n x 1)
%   Sigma0  initial prior covariance (n x n)
%   Z       pruning parameter (i.e. number of most likely Gaussian 
%           components to keep at each step)
%
% OUTPUTS:
%   mu      posterior mean (n x T matrix of n x 1 vectors)
%   Sigma   posterior covariance (n x n x T array of n x n matrices)
%   alpha   mixture weights (T x 1)
%
%==========================================================================
function [mu,Sigma,alpha] = MHKF(A,B,C,Q,R,u,y,mu0,Sigma0,beta,gamma,M,L,Z)
    
    % length of time vector
    T = length(u);
    
    % number of Gaussian components in initial prior
    N0 = length(mu0);
    
    % determines n (dimension of state vector)
    n = length(mu0{1});
    
    % identity matrix for update step
    I = eye(n);
    
    % preallocates arrays
    mu = zeros(n,T,Z);
    Sigma = zeros(n,n,T,Z);
    alpha = zeros(Z,T);
    
    % assigns initial conditions
    for i = 1:N0
        mu(:,1,i) = mu0{i};
        Sigma(:,:,1,i) = Sigma0{i};
        alpha(i,1) = 1/N0;
    end
    
    % array to store number of Gaussian components in prior at each t
    N = zeros(T,1);
    N(1) = N0;
    
    % filtering
    for t = 2:T
        
        % preallocates arrays for predict step
        mu_predict = zeros(n,N(t-1),M);
        Sigma_predict = zeros(n,n,N(t-1),M);
        alpha_predict = zeros(N(t-1),M);
        
        % predict step
        for i = 1:N(t-1)
            for j = 1:M
                
                % calculates necessary quantities
                mui = mu(:,t-1,i);
                Aj = A(mui,u(:,t-1),t,j);
                Bj = B(mui,u(:,t-1),t,j);
                
                % predicts mean, covariance, and corresponding weight
                mu_predict(:,i,j) = Aj*mui+Bj*u(:,t-1);
                Sigma_predict(:,:,i,j) = Aj*Sigma(:,:,t-1,j)*Aj'+Q(j);
                alpha_predict(i,j) = beta(j)*alpha(i,t-1);
                
            end
        end
        
        % reindex results of predict step
        mu_predict = mu_predict(:,:);
        Sigma_predict = Sigma_predict(:,:,:);
        alpha_predict = alpha_predict(:);
        
        % preallocates arrays for update step
        mu_update = zeros(n,L,M*N(t-1));
        Sigma_update = zeros(n,n,L,M*N(t-1));
        alpha_bar = zeros(L,M*N(t-1));
        
        % update step
        for i = 1:L
            for k = 1:(M*N(t-1))
                
                % calculates/extracts required quantities
                muk = mu_predict(:,k);
                Sigmak = Sigma_predict(:,:,k);
                yt = y(:,t);
                Ci = C(muk,t,i);
                
                % Kalman gain
                K = Sigmak*Ci'*inv(Ci*Sigmak*Ci'+R(i));
                
                % updates mean and covariance
                mu_update(:,i,k) = muk+K*(yt-Ci*muk);
                Sigma_update(:,:,i,k) = (I-K*Ci)*Sigmak;
                
                % updates unnormalized weights
                muy = Ci*muk;
                Sigmay = Ci*Sigmak*Ci'+R(i);
                pik = mvnpdf(yt,muy,Sigmay);
                alpha_bar(i,k) = gamma(i)*alpha_predict(k)*pik;
                
            end
        end
        
        % normalizes weights
        S = 0;
        for i = 1:L
            for k = 1:(M*N(t-1))
                S = S+alpha_bar(i,k);
            end
        end
        alpha_update = alpha_bar/S;
        
        % reindex results of update step
        mu_update = mu_update(:,:);
        Sigma_update = Sigma_update(:,:,:);
        alpha_update = alpha_update(:);
        
        % finds indices of the Z most likely Gaussian components
        [~,iZ] = maxk(alpha_update,Z);
        
        % prunes results of update step
        mu_update = mu_update(:,iZ);
        Sigma_update = Sigma_update(:,:,iZ);
        alpha_update = alpha_update(iZ);
        
        % stores results
        mu(:,t,:) = mu_update;
        Sigma(:,:,t,:) = Sigma_update;
        alpha(:,t) = alpha_update;
        
        % number of Gaussian components
        N(t) = size(mu_update,2);
        
    end
    
end