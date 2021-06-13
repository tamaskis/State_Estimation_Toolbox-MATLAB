%==========================================================================
%
% FastSLAM  FastSLAM algorithm for simultaneous localization and mapping.
%
%   [mu,Sigma] = FastSLAM(f,g,w,pv,u,y,x0,N)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-06-03
%
%==========================================================================
%
% INPUTS:
%   w_noise         (3 x t) time history of process noise
%   R               (2 x 2) measurement noise covariance
%   u               (m x T matrix of m x 1 vectors) control input
%   y               (k x (T-1) matrix of k x 1 vectors) measurement
%   mu0_pose        (3 x 1) initial prior mean for pose
%   mu0_feature     (2 x M) initial prior mean for features
%   Sigma0_pose     (3 x 3) initial prior covariance for pose
%   Sigma0_feature  (3 x 3) initial prior covariance for features
%   N               (1 x 1) number of particles
%   dt              (1 x 1) time step [s]
%   i_closest       (10 x 1) indices of the 10 closest features
%
% OUTPUTS:
%   mu_pose         (3 x T) posterior pose means
%   Sigma_pose      (n x n x T) posterior pose covariances
%   mu_feature      (2 x M x N x T) posterior feature location means
%   Sigma_feature   (2 x 2 x M x N x T) posterior feature location
%                                       covariances
%
%==========================================================================
function [mu_pose,Sigma_pose,mu_feature,Sigma_feature] =...
    FastSLAM(w_noise,R,u,y,mu0_pose,mu0_feature,Sigma0_feature,M,N,dt,...
    i_closest)
    
    % determines length of time vector
    T = length(u);
    
    % preallocates arrays for pose and pose statistics
    mu_pose = zeros(3,T);           % pose estimate means
    Sigma_pose = zeros(3,3,T);      % pose estimate covariances
    x_predict = zeros(3,N);      	% robot pose
    
    % preallocates arrays for feature statistics
    mu_update = zeros(2,M,N);
    mu_feature = zeros(2,M,N,T);
    Sigma_update = zeros(2,2,M,N);
    Sigma_feature = zeros(2,2,M,N,T);
    
    % initial particle set
    x = zeros(3,1,N);
    for i = 1:N
        x(:,i) = mu0_pose;
        mu_feature(:,:,i,1) = mu0_feature;
        Sigma_feature(:,:,:,i,1) = Sigma0_feature;
    end
    
    % initial weights
    w = 1/N*ones(N,1);
    
    % filtering
    for t = 2:T
        
        % unpacks control vector
        nu = u(1,t-1);
        phi = u(2,t-1);
            
        % prediction step
        for i = 1:N
            
            % extracts pose from state vector
            px = x(1,i);
            py = x(2,i);
            theta = x(3,i);
            
            % propagates pose
            x_predict(:,i) = [px+nu*cos(theta)*dt;py+nu*sin(theta)*dt;...
                theta+phi*dt]+w_noise(:,t-1);
        end
        
        % preallocates w_hat
        w_hat = zeros(N,1);
        
        % update step
        for i = 1:N
            
            % initializes w_bar
            w_bar = zeros(M,1);
            
            for j = 1:M
                
                % extracts mean and covariance of specific feature at
                % previous time step
                mu_ij = mu_feature(:,j,i,t-1);
                Sigma_ij = Sigma_feature(:,:,j,i,t-1);

                % determines C matrix, true measurement, and expected
                % measurement
                if ismember(j,i_closest(t))
                    
                    % x and y coordinates of jth feature location for ith
                    % particle
                    mx = mu_ij(1);
                    my = mu_ij(2);
                    
                    % robot pose from prediction step for ith particle
                    px = x_predict(1,i);
                    py = x_predict(2,i);
                    theta = x_predict(3,i);
                    
                    % C matrix
                    C = [mx/sqrt((mx-px)^2+(my-py)^2)   my/sqrt((mx-px)^2+(my-py)^2);
                        (py-my)/((mx-px)^2+(my-py)^2)   (mx-px)/((mx-px)^2+(my-py)^2)];
                    
                    % finds the index of the measurement corresponding to
                    % the closest feature
                    k = find(i_closest(t)==j);
                    
                    % true and predicted measurements
                    yj = y((2*k-1):(2*k),t);
                    y_exp = [sqrt((mx-px)^2+(my-py)^2);atan2(mx-py,my-...
                        px)-theta];
                    
                else
                    C = zeros(2,2);
                    yj = [0;0];
                    y_exp = [0;0];
                end
                
                % Kalman gain
                K = Sigma_ij*C'*inv(C*Sigma_ij*C'+R);
                
                % residual (wrapping to pi where needed)
                z = yj-y_exp;
                z(2) = wrapToPi(z(2));
                
                % updates feature mean and covariance
                mu_update(:,j,i) = mu_ij+K*z;
                Sigma_update(:,:,j,i) = (eye(2)-K*C)*Sigma_ij;
                
                % calculates weight
                w_bar(j) = mvnpdf(yj,y_exp,inv(C*Sigma_ij*C'+R));
                
            end
            
            % un-normalized weight
            w_hat(i) = prod(w_bar)*w(i);
            
        end
        
        % updated weights
        w = w_hat/sum(w_hat);
        
        % resample step
        for j = 1:N
            s = rand;
            w_cumulative_sum = 0;
            i = 1;
            new_particle = false;
            while (i < N) && (new_particle == false)
                w_cumulative_sum = w_cumulative_sum+w(i);
                if s < w_cumulative_sum
                    x(:,j) = x_predict(:,i);
                    mu_feature(:,:,j,t) = mu_update(:,:,i);
                    Sigma_feature(:,:,:,j,t) = Sigma_update(:,:,:,i);
                    new_particle = true;
                end
                i = i+1;
            end
        end
        
        % resets the weights to 1/N
        w(:) = 1/N;
        
        % posterior (empirical) mean for pose
        for i = 1:N
            mu_pose(:,t) = mu_pose(:,t)+w(i)*x(:,i);
        end
        
        % posterior (empirical) covariance for pose
        for i = 1:N
            Sigma_pose(:,:,t) = Sigma_pose(:,:,t)+w(i)*(x(:,i)-...
                mu_pose(:,t))*(x(:,i)-mu_pose(:,t))';
        end

    end
    
end