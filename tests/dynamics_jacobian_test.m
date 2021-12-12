%% dynamics_jacobian_test.m
% Unit testing of the dynamics_jacobian function.
%
% Author: Tamas Kis
% Last Update: 2021-11-30



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;

% adds path to all "State Estimation" functions
addpath(genpath('../src'));
addpath(genpath('../res'));



%% NOTE

% This example is taken from AA 273 Problem Set 4 (Stanford University).



%% TEST

% time step [s]
dt = 0.001;

% principal inertia tensor [kg.m^2]
J = [1  0  0;
     0  5  0;
     0  0  5];

% discrete-time nonlinear dynamics equation
fd = @(x,u,k) discrete_dynamics(x,u,J,dt);

% test condition
x = [10;11;12];
u = [3;2;1];
k = 20;

% evaluates dynamics Jacobian numerically
actual = dynamics_jacobian(fd,x,u,k);
expected = exact_dynamics_jacobian(x,J,dt);

% unit test
TEST_EQUAL(actual,expected);


% continuous-time nonlinear dynamics equation
f = @(x,u,k) continuous_dynamics(x,u,J);
actual = eye(3)+dynamics_jacobian(f,x,u,k)*dt;

% unit test
TEST_EQUAL(actual,expected);



%% ADDITIONAL FUNCTIONS

%==========================================================================
% continuous_dynamics  Continuous-time nonlinear dynamics equation.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   omega       - (3×1 double) angular velocity (state vector) [rad/s]
%   tau         - (3×1 double) control torque (control input) [N.m]
%   J           - (3×3 double) principal inertia tensor [kg.m^2]
%
% -------
% OUTPUT:
% -------
%   domega      - (3×1 double) angular acceleration (state vector 
%                 derivative) [rad/s^2]
%
%==========================================================================
function domega = continuous_dynamics(omega,tau,J)

    % unpacks angular velocity vector
    wx = omega(1);
    wy = omega(2);
    wz = omega(3);

    % unpacks control input
    taux = tau(1);
    tauy = tau(2);
    tauz = tau(3);
    
    % unpacks principal inertia tensor
    Jx = J(1,1);
    Jy = J(2,2);
    Jz = J(3,3);
    
    % evalutes f(x,u,t)
    domega = [((Jy-Jz)*wy*wz+taux)/Jx;
              ((Jz-Jx)*wz*wx+tauy)/Jy;
              ((Jx-Jy)*wx*wy+tauz)/Jz];

end




%==========================================================================
% discrete_dynamics  Discrete-time nonlinear dynamics equation.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   omega       - (3×1 double) angular velocity (state vector) at current
%                 sample time [rad/s]
%   tau         - (3×1 double) control torque (control input) at current 
%                 sample time [N.m]
%   J           - (3×3 double) principal inertia tensor [kg.m^2]
%   dt          - (1×1 double) time step [s]
%
% -------
% OUTPUT:
% -------
%   omega_new   - (3×1 double) angular velocity (state vector) at next
%                 sample time [rad/s]
%
%==========================================================================
function omega_new = discrete_dynamics(omega,tau,J,dt)

    % unpacks angular velocity vector
    wx = omega(1);
    wy = omega(2);
    wz = omega(3);
    
    % unpacks principal inertia tensor
    Jx = J(1,1);
    Jy = J(2,2);
    Jz = J(3,3);
    
    % evalutes f(x,u,t)
    omega_new = [wx+((Jy-Jz)*wy*wz*dt)/Jx;
                 wy+((Jz-Jx)*wz*wx*dt)/Jy;
                 wz+((Jx-Jy)*wx*wy*dt)/Jz]+[dt/Jx,dt/Jy,dt/Jz]*tau;
    
end



%==========================================================================
% exact_dynamics_jacobian  Exact dynamics Jacobian.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   omega   - (3×1 double) angular velocity (state vector) at current
%             sample time [rad/s]
%   J       - (3×3 double) principal inertia tensor [kg.m^2]
%   dt      - (1×1 double) time step [s]
%
% -------
% OUTPUT:
% -------
%   F       - (3×3 double) dynamics Jacobian at current sample time
%
%==========================================================================
function F = exact_dynamics_jacobian(omega,J,dt)

    % unpacks angular velocity vector [rad/s]
    wx = omega(1);
    wy = omega(2);
    wz = omega(3);
    
    % unpacks principal inertia tensor [kg.m^2]
    Jx = J(1,1);
    Jy = J(2,2);
    Jz = J(3,3);
    
    % assembles dynamics Jacobian
    F = [1                  (Jy-Jz)*wz*dt/Jx   (Jy-Jz)*wy*dt/Jx;
         (Jz-Jx)*wz*dt/Jy   1                  (Jz-Jx)*wx*dt/Jy;
         (Jx-Jy)*wy*dt/Jz   (Jx-Jy)*wx*dt/Jz   1];
    
end