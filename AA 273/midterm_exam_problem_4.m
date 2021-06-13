%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Midterm Exam Problem 4(c)

%-------------------------------------------------------------------------%



%% SCRIPT SETUP

% clears variables and command window, closes all figures
clear;
clc;
close all;

% plot parameters
plot_position = [540,300,700,500]; % plot position [x,y,l,w]
line_width = 1.5; % line width [#]
axis_font_size = 18; % axis label font size [#]
legend_font_size = 14; % legend font size [#]
title_font_size = 20; % title font size [#]



%% TRUE TRAJECTORY AND MEASUREMENT PROCESS SIMULATIONS

% initial condition
x0 = [0;70;10];

% time parameters
dt = 1; % time step [s]
tf = 100*dt; % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% system matrices
A = [1,dt,0;0,1,0;0,0,1];
C = [1,0,0;0,1,1];

% dimension parameters
n = size(A,1);
k = size(C,1);

% noise covariances
q = 3;
Q = [0,0,0;0,q,0;0,0,0];
R = [5,0;0,5];

% simulated noise (with 0 mean)
w = gaussian_random_sample(zeros(1,1),q,T-1);
v = gaussian_random_sample(zeros(2,1),R,T-1);

% array preallocation
x = zeros(n,T);
y = zeros(k,T);

% stores initial condition
x(:,1) = x0;

% simulation
for tt = 1:(T-1)
    x(:,tt+1) = A*x(:,tt)+[0;w(:,tt);0];
    y(:,tt+1) = C*x(:,tt+1)+v(:,tt);
end

% extracts true position, velocity, and bias from state vector
P_true = x(1,:);
S_true = x(2,:);
dS_true = x(3,:);



%% KALMAN FILTERING

% initial state estimate
mu0 = [0;60;0];
Sigma0 = [1,0,0;0,10,0;0,0,100];

% identity matrix for update step
I = eye(n);

% preallocates arrays
mu = zeros(n,T);
Sigma = zeros(n,n,T);

% assigns initial conditions
mu(:,1) = mu0;
Sigma(:,:,1) = Sigma0;

% runs Kalman filter until tf
for tt = 2:T

    % prediction step
    mu_predict = A*mu(:,tt-1);
    Sigma_predict = A*Sigma(:,:,tt-1)*A'+Q;

    % Kalman gain
    K = Sigma_predict*C'*inv(C*Sigma_predict*C'+R);

    % update step
    mu(:,tt) = mu_predict+K*(y(:,tt)-C*mu_predict);
    Sigma(:,:,tt) = (I-K*C)*Sigma_predict;

end

% extracts estimated position, velocity, and bias from estimated state
P_est = mu(1,:);
S_est = mu(2,:);
dS_est = mu(3,:);



%% PLOTS

% true and estimated position
figure('position',plot_position);
hold on;
plot(t,P_est,'linewidth',line_width);
plot(t,P_true,'linewidth',line_width);
hold off;
grid on;
xlabel('$t$','interpreter','latex','fontsize',axis_font_size);
ylabel('$P_{t}$','interpreter','latex','fontsize',axis_font_size);
title('\textbf{Estimated and True Position}','interpreter','latex',...
    'fontsize',title_font_size);
legend('estimated position','true position','interpreter','latex',...
    'fontsize',legend_font_size,'location','best');

% true and estimated velocity
figure('position',plot_position);
hold on;
plot(t,S_est,'linewidth',line_width);
plot(t,S_true,'linewidth',line_width);
hold off;
grid on;
xlabel('$t$','interpreter','latex','fontsize',axis_font_size);
ylabel('$S_{t}$','interpreter','latex','fontsize',axis_font_size);
title('\textbf{Estimated and True Velocity}','interpreter','latex',...
    'fontsize',title_font_size);
legend('estimated velocity','true velocity','interpreter','latex',...
    'fontsize',legend_font_size,'location','best');

% true and estimated bias
figure('position',plot_position);
hold on;
plot(t,dS_est,'linewidth',line_width);
plot(t,dS_true,'linewidth',line_width);
hold off;
grid on;
xlabel('$t$','interpreter','latex','fontsize',axis_font_size);
ylabel('$\Delta S$','interpreter','latex','fontsize',axis_font_size);
title('\textbf{Estimated and True Bias}','interpreter','latex',...
    'fontsize',title_font_size);
legend('estimated bias','true bias','interpreter','latex','fontsize',...
    legend_font_size,'location','best');



%% MEAN AND VARIANCE OF YT-ST

% extracts mu_Delta and sigma_Delta from most recent estimate
mu_Delta = mu(3,end)
sigma_Delta2 = Sigma(3,3,end)

% mean and variance of yt-st
mu_ytst = mu_Delta
sigma_ytst2 = sigma_Delta2+5


%% PROBABILITY CALCULATIONS

ytst = norminv(0.05,10,3)



%% ADDITIONAL FUNCTIONS

% INPUT: mu - mean (n-by-1 vector)
%        Sigma - covariance matrix (n-by-n matrix)
%        N - number of samples
% OUTPUT: x - sample points (n-by-N matrix, each column stores one vector
function x = gaussian_random_sample(mu,Sigma,N)
    
    % determines dimension of random vector
    n = length(mu);
    
    % random samples from standard normal distribution
    w = randn(N,n);
    
    % converts to Gaussian distribution with given parameters
    x = (w*chol(Sigma)+repmat(mu',N,1))';
    
end