%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Problem Set 7 Problem 2

%-------------------------------------------------------------------------%



%% SCRIPT SETUP

% clears variables and command window, closes all figures
clear;
clc;
close all;

% adds path to Estimation Toolbox
addpath("..");

% loads plot parameters
PLOT_PARAMETERS;

% seeds random number generators
rng(2);



%% MAP

% coordinates of upper half
x_upper = [1.91;1.17;-0.42;-0.42;0.42;0.42;1.91;1.91;1.07;-1.07;-1.91;...
    -1.91];
y_upper = [0;0.7;0.7;1.49;1.49;0.96;0.96;2.11;2.93;2.93;2.11;0];

% coordinates of lower half
x_lower = -x_upper;
y_lower = -y_upper;

% concatenates coordinates to create "outer" line
x_outer = [x_upper;x_lower];
y_outer = [y_upper;y_lower];

% offsets for upper half of "inner" line
o = 0.1;
x_offsets = [-o;-o*cosd(60);-o;-o;o;o;-o;-o;-o*cosd(60);o*cosd(60);o;o];
y_offsets = [-o*sind(30);-o;-o;o;o;o;o;-o*sind(30);-o;-o;-o*sind(30);o*...
    sind(30)];

% offsets for both halves of "inner" line
x_offsets = [x_offsets;-x_offsets];
y_offsets = [y_offsets;-y_offsets];

% coordinates for "inner" line
x_inner = x_outer+x_offsets;
y_inner = y_outer+y_offsets;

% increases number of points through interpolation
[x_outer,y_outer] = interpolate_shape(x_outer,y_outer,1);
[x_inner,y_inner] = interpolate_shape(x_inner,y_inner,1);

% plots Stanford S
figure('position',plot_position);
hold on;
plot(x_outer,y_outer,'k','linewidth',3);
plot(x_inner,y_inner,'k','linewidth',3);
hold off;
axis equal;

% groups all feature points together and adds six extra points
x = [x_outer;x_inner;-0.05;0.05;0.05;-0.05;-0.05;0.05];
y = [y_outer;y_inner;-0.1;-0.1;0.1;0.1;0;0];
figure('position',plot_position);
scatter(x,y);axis equal;

% creates (2 x 100) matrix of true feature locations
m_true = [x';y'];

% scales the feature locations
m_true = m_true*10;
x_outer = x_outer*10;
x_inner = x_inner*10;
y_outer = y_outer*10;
y_inner = y_inner*10;



%% SIMULATION OF NONLINEAR DISCRETE TIME SYSTEM

% time parameters
dt = 0.1; % time step [s]
tf = 50; % simulation end time [s]

% time vector and its length
t = 0:dt:tf;
T = length(t);

% dimension parameters
n = 203;
k = 20;

% noise covariances
Q = [0.1*eye(3)*dt,zeros(3,200);zeros(200,3),0.0001*eye(200)];
R = 0.1*eye(k);

% initial state estimate (prior distribution)
mu0 = [0;0;0;m_true(:)+0.2];
Sigma0 = 0.01*eye(length(mu0));

% simulated noise
w = gaussian_random_sample(zeros(n,1),Q,T);
v = gaussian_random_sample(zeros(k,1),R,T);

% sets process noise samples of feature locations to 0
w(4:203,:) = 0;

% control signal
u = [ones(size(t));sin(t)];

% array preallocation
x = zeros(n,T);
y = zeros(k,T);
i_closest = zeros(10,T);

% samples initial condition from the prior distribution
x(:,1) = gaussian_random_sample(mu0,Sigma0);

% initial closest features
i_closest(:,1) = closest_features(x(1:2,1),m_true);

% functions for nonlinear dynamics and nonlinear measurements
f = @(x,u,t) f_discrete(x,u,t,dt);
g = @(x,t,i_closest) g_discrete(x,t,i_closest);

% simulation (using true feature locations)
for tt = 1:(T-1)
    
    % propagates state vector with noise
    x(:,tt+1) = f([x(1:3,tt);m_true(:)],u(:,tt),tt)+w(:,tt);
    
    % indices of closest features
    i_closest(:,tt+1) = closest_features(x(1:2,tt+1),m_true);
    
    % measurement with noise
    y(:,tt+1) = g([x(1:3,tt+1);m_true(:)],tt+1,i_closest)+...
        v(:,tt+1);
    
end

% extracts time history of true position
px_true = x(1,:);
py_true = x(2,:);

% time history of true closest feature locations
mx_true = zeros(10,T);
my_true = zeros(10,T);
for i = 2:T
    m = x(4:203,i);
    m = reshape(m,2,100);
    m_closest = m(:,i_closest(:,i));
    mx_true(:,i) = m_closest(1,:)';
    my_true(:,i) = m_closest(2,:)';
end

% plots Stanford S with feature locations
figure('position',plot_position);
hold on;
plot(px_true,py_true,'linewidth',line_width,'color',cardinal_red);
plot(x_outer,y_outer,'k:','linewidth',1.5);
plot(x_inner,y_inner,'k:','linewidth',1.5);
scatter(m_true(1,:),m_true(2,:),25,'k','filled');
scatter(mx_true,my_true,50,matlab_blue,'linewidth',line_width);
hold off;
axis equal;



%% EKF SLAM

% dynamics and measurement Jacobians
A = @(x,u,t) dynamics_jacobian(x,u,t,dt);
C = @(x,t) measurement_jacobian(x,t,i_closest);

% runs EKF SLAM algorithm
tic
[mu,Sigma] = EKF_SLAM(f,g,A,C,Q,R,u,y,mu0,Sigma0,i_closest);
toc

% extracts time history of position estimates
px_ekfslam = mu(1,2:end);
py_ekfslam = mu(2,2:end);

% extracts time history of closest feature location estimates
mx_ekfslam = zeros(10,T);
my_ekfslam = zeros(10,T);
for i = 2:T
    m = mu(4:203,i);
    m = reshape(m,2,100);
    m_closest = m(:,i_closest(:,i));
    mx_ekfslam(:,i) = m_closest(1,:)';
    my_ekfslam(:,i) = m_closest(2,:)';
end

% initializes figure
figure('position',two_subplot_position);

% true trajectory and map
subplot(1,2,1);
hold on;
plot(x_outer,y_outer,'k:','linewidth',1.5,'handlevisibility','off');
plot(x_inner,y_inner,'k:','linewidth',1.5,'handlevisibility','off');
scatter(m_true(1,:),m_true(2,:),25,'k','filled');
axis equal;
plot(px_true,py_true,'linewidth',line_width,'color',matlab_blue);
xlabel('$x$','interpreter','latex','fontsize',axis_font_size);
ylabel('$y$','interpreter','latex','fontsize',axis_font_size);
title('\textbf{True Trajectory and Map}','interpreter','latex',...
    'fontsize',title_font_size);
legend('true feature locations','true trajectory','interpreter','latex',...
    'fontsize',legend_font_size);
hold off;

% EKF SLAM results
subplot(1,2,2);
hold on;
plot(x_outer,y_outer,'k:','linewidth',1.5,'handlevisibility','off');
plot(x_inner,y_inner,'k:','linewidth',1.5,'handlevisibility','off');
axis equal;
scatter(mx_ekfslam(:),my_ekfslam(:),60,'g','linewidth',line_width);
plot(px_ekfslam,py_ekfslam,'linewidth',line_width,'color',cardinal_red);
xlabel('$x$','interpreter','latex','fontsize',axis_font_size);
ylabel('$y$','interpreter','latex','fontsize',axis_font_size);
title('\textbf{EKF SLAM Results}','interpreter','latex','fontsize',...
    title_font_size);
legend('estimated feature locations','estimated trajectory',...
    'interpreter','latex','fontsize',legend_font_size);
hold off;



%% FAST SLAM

% number of particles to use
N = 10;

% number of features
M = 100;

% measurement noise for feature locations
R = 0.1*eye(2);

% initial prior distribution for pose
mu0_pose = mu0(1:3);
Sigma0_pose = Sigma0(1:3,1:3);

% initial prior distribution for feature locations
mu0_feature = m_true+0.2;
Sigma0_feature = zeros(2,2,M);
for i = 1:M
    Sigma0_feature(:,:,i) = 0.01*eye(2);
end

% performs FastSLAM
tic
[mu_pose,Sigma_pose,mu_feature] = FastSLAM(w(1:3,:),R,u,y,mu0_pose,...
    mu0_feature,Sigma0_feature,M,N,dt,i_closest);
toc

% extracts time history of position estimates
px_fastslam = mu_pose(1,2:end);
py_fastslam = mu_pose(2,2:end);

% extracts time history of closest feature location estimates
mx_fastslam = zeros(10,T,N);
my_fastslam = zeros(10,T,N);
for i = 2:T
    for j = 1:N
        m_closest = mu_feature(:,i_closest(:,i),j,i);
        mx_fastslam(:,i,j) = m_closest(1,:)';
        my_fastslam(:,i,j) = m_closest(2,:)';
    end
end

% extracts just the first particle
mx_fastslam = mx_fastslam(:,:,1);
my_fastslam = my_fastslam(:,:,1);

% initializes figure
figure('position',two_subplot_position);

% true trajectory and map
subplot(1,2,1);
hold on;
plot(x_outer,y_outer,'k:','linewidth',1.5,'handlevisibility','off');
plot(x_inner,y_inner,'k:','linewidth',1.5,'handlevisibility','off');
scatter(m_true(1,:),m_true(2,:),25,'k','filled');
axis equal;
plot(px_true,py_true,'linewidth',line_width,'color',matlab_blue);
xlabel('$x$','interpreter','latex','fontsize',axis_font_size);
ylabel('$y$','interpreter','latex','fontsize',axis_font_size);
title('\textbf{True Trajectory and Map}','interpreter','latex',...
    'fontsize',title_font_size);
legend('true feature locations','true trajectory','interpreter','latex',...
    'fontsize',legend_font_size);
hold off;

% FastSLAM results
subplot(1,2,2);
hold on;
plot(x_outer,y_outer,'k:','linewidth',1.5,'handlevisibility','off');
plot(x_inner,y_inner,'k:','linewidth',1.5,'handlevisibility','off');
axis equal;
scatter(mx_fastslam(:),my_fastslam(:),60,'g','linewidth',line_width);
plot(px_fastslam,py_fastslam,'linewidth',line_width,'color',cardinal_red);
xlabel('$x$','interpreter','latex','fontsize',axis_font_size);
ylabel('$y$','interpreter','latex','fontsize',axis_font_size);
title('\textbf{FastSLAM Results}','interpreter','latex','fontsize',...
    title_font_size);
legend('estimated feature locations','estimated trajectory',...
    'interpreter','latex','fontsize',legend_font_size);
hold off;



%% ADDITIONAL FUNCTIONS

%==========================================================================
% Closest features.
%==========================================================================
%
% INPUTS:
%   p       (2 x 1) position
%   m_true  (2 x 100) true feature locationss
%
% OUTPUTS:
%   i_closest   (10 x 1) indices of the 10 closest features
%
%==========================================================================
function i_closest = closest_features(p,m_true)
    
    % distances between positions and feature locations
    s = zeros(100,1);
    for i = 1:100
        s(i) = norm(m_true(:,i)-p);
    end
    
    % adds a column of indices
    s = [(1:100)',s];
    
    % reorders s from smallest to largest position
    s = sortrows(s,2);
    
    % indices of the 10 closest features
    i_closest = s(1:10,1);
    
end



%==========================================================================
% Discrete-time nonlinear dynamics.
%==========================================================================
%
% INPUTS:
%   x       (203 x 1) state vector at current time step
%   u       (2 x 1) control input at current time step
%   t       (1 x 1) current iteration (corresponding to discrete time)
%   dt      (1 x 1) [s] time step
%
% OUTPUTS:
%   f_eval  (203 x 1) evaluation of f (produces state vector at next time 
%                     step)
%
%==========================================================================
function f_eval = f_discrete(x,u,t,dt)

    % extracts pose from state vector
    px = x(1);
    py = x(2);
    theta = x(3);
    
    % unpacks control vector
    nu = u(1);
    phi = u(2);
    
    % evalutes f(xt,ut)
    f_eval = [px+nu*cos(theta)*dt;
              py+nu*sin(theta)*dt;
              theta+phi*dt;
              x(4:203)];
    
end



%==========================================================================
% Discrete-time nonlinear measurement.
%==========================================================================
%
% INPUTS:
%   x           (203 x 1) state vector at current time step
%   t           current iteration (corresponding to discrete time)
%   i_closest   (10 x 1) indices of 10 closest features
%
% OUTPUTS:
%   g_eval      evaluation of g (produces measurement at current time step)
%
%==========================================================================
function g_eval = g_discrete(x,t,i_closest)
    
    % extracts pose from state vector
    px = x(1);
    py = x(2);
    theta = x(3);
    
    % extracts feature locations from state vector
    m = x(4:203);
    
    % reshapes m to be (2 x 100)
    m = reshape(m,2,100);
    
    % extracts estimates of the 10 closest feature locations
    m = m(:,i_closest(:,t));
    
    % evalutes g(xt)
    g_eval = zeros(20,1);
    for i = 1:10
        g_eval(i) = sqrt((m(1,i)-px)^2+(m(2,i)-py)^2);
    end
    for i = 11:20
        g_eval(i) = atan2(m(2,i-10)-py,m(1,i-10)-px)-theta;
    end
                      
end



%==========================================================================
% Dynamics Jacobian.
%==========================================================================
%
% INPUTS:
%   x       (203 x 1) state vector
%   u       (2 x 1) control input
%   t       (1 x 1) current iteration (corresponding to discrete time)
%   dt      (1 x 1) [s] time step
%
% OUTPUTS:
%   A       (203 x 203) dynamics Jacobian
%
%==========================================================================
function A = dynamics_jacobian(x,u,t,dt)

    % extracts "theta" from state vector
    theta = x(3);
    
    % extracts "nu" from control input
    nu = u(1);
    
    % A1
    A1 = [1   0   -nu*sin(theta)*dt;
          0   1    nu*cos(theta)*dt;
          0   0    1];
      
    % A2
    A2 = eye(200);
    
    % assembles full dynamics Jacobian
    A = [A1,zeros(3,200);zeros(200,3),A2];
    
end



%==========================================================================
% Measurement Jacobian.
%==========================================================================
%
% INPUTS:
%   x       (203 x 1) state vector
%   t       (1 x 1) current iteration (corresponding to discrete time)
%   i_closest   (10 x 1) indices of 10 closest features
%
% OUTPUTS:
%   C       (20 x 203) measurement Jacobian
%
%==========================================================================
function C = measurement_jacobian(x,t,i_closest)
    
    % extracts position from state vector
    px = x(1);
    py = x(2);
    
    % extracts feature locations from state vector
    m = x(4:203);
    
    % reshapes m to be (2 x 100)
    m = reshape(m,2,100);
    
    % extracts estimates of the 10 closest feature locations
    m = m(:,i_closest(:,t));
    
    % C1 matrix
    C1 = zeros(10,3);
    for i = 1:10
        mx = m(1,i);
        my = m(2,i);
        C1(i,:) = [(px-mx)/sqrt((mx-px)^2+(my-py)^2),(py-my)/sqrt((mx-...
            px)^2+(my-py)^2),0];
    end
    
    % C2 matrix
    C2 = zeros(10,200);
    for i = 1:10
        for j = 1:100
            if (i_closest(i) == j)
                mx = m(1,i);
                my = m(2,i);
                C2(i,2*j-1) = mx/sqrt((mx-px)^2+(my-py)^2);
                C2(i,2*j) = my/sqrt((mx-px)^2+(my-py)^2);
            else
                C2(i,2*j-1) = 0;
                C2(i,2*j) = 0;
            end
        end
    end
    
    % C3 matrix
    C3 = zeros(10,3);
    for i = 1:10
        mx = m(1,i);
        my = m(2,i);
        C3(i,:) = [(my-py)/sqrt((mx-px)^2+(my-py)^2),(px-mx)/sqrt((mx-...
            px)^2+(my-py)^2),-1];
    end
    
    % C4 matrix
    C4 = zeros(10,200);
    for i = 1:10
        for j = 1:100
            if (i_closest(i) == j)
                mx = m(1,i);
                my = m(2,i);
                C4(i,2*j-1) = (py-my)/((mx-px)^2+(my-py)^2);
                C4(i,2*j) = (mx-px)/((mx-px)^2+(my-py)^2);
            else
                C4(i,2*j-1) = 0;
                C4(i,2*j) = 0;
            end
        end
    end
    
    % assembles full measurement Jacobian
    C = [C1   C2;
         C3   C4];
    
end




















%==========================================================================
% State transition matrix.
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   x       - (4×1) state vector
%   u       - (2×1) control input
%   t    	- (1×1) current iteration (corresponding to discrete time)
%   dt      - (1×1) time step
%
% --------
% OUTPUTS:
% --------
%   A       - (4×4) state transition matrix
%
%==========================================================================
function A = state_transition_matrix(x,u,t,dt)
    A = [eye(2),eye(2)*dt;zeros(2),eye(2)];
end



%==========================================================================
% Input matrix.
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   x       - (4×1) state vector
%   u       - (2×1) control input
%   t    	- (1×1) current iteration (corresponding to discrete time)
%   dt      - (1×1) time step
%
% --------
% OUTPUTS:
% --------
%   B       - (4×2) input matrix
%
%==========================================================================
function B = input_matrix(x,u,t,dt)
    B = [zeros(2);eye(2)*dt];
end



%==========================================================================
% Observation matrix.
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   x       - (4×1) state vector
%   u       - (2×1) control input
%   t    	- (1×1) current iteration (corresponding to discrete time)
%
% --------
% OUTPUTS:
% --------
%   C   	- (4×2) observation matrix
%
%==========================================================================
function C = observation_matrix(x,u,t)
    C = [eye(2),zeros(2)];
end