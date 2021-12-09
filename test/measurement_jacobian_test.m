%% measurement_jacobian_test.m
% Unit testing of the measurement_jacobian function.
%
% Author: Tamas Kis
% Last Update: 2021-11-13



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

% test condition
x = [5;6;7];
u = [3;2;1];
k = 20;
c = 10;

% principal inertia tensor [kg.m^2]
J = [1  0  0;
     0  5  0;
     0  0  5];

% discrete-time nonlinear measurement equation
hd = @(x,k) sat(x,c);

% evaluates dynamics Jacobian numerically
obtained = measurement_jacobian(hd,x,k);
expected = exact_measurement_jacobian(x,c);

% unit test
TEST_EQUALITY(obtained,expected);



%% ADDITIONAL FUNCTIONS

%==========================================================================
% satx  Element-wise saturation function (discrete-time nonlinear 
% measurement equation).
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (n×1 double) independent variable
%   c       - (1×1 double) threshold value for saturation
%
% -------
% OUTPUT:
% -------
%   satx    - (n×1 double) evaluation of sat(x)
%
% -----
% NOTE:
% -----
%   --> n = dimension of independent variable
%
%==========================================================================
function satx = sat(x,c)
    satx = zeros(size(x));
    for i = 1:length(x)
        if iabs(x(i)) < c
            satx(i) = x(i);
        else
            satx(i) = c*sign(x(i));
        end
    end
end



%==========================================================================
% exact_measurement_jacobian  Exact measurement Jacobian.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   omega   - (3×1 double) angular velocity (state vector) at current 
%             sample time [rad/s]
%   c       - (1×1 double) threshold value for measurement saturation 
%             [rad/s]
%
% -------
% OUTPUT:
% -------
%   H       - (3×3 double) measurement Jacobian at current sample time
%
%==========================================================================
function C = exact_measurement_jacobian(omega,c)

    % unpacks angular velocity
    wx = omega(1);
    wy = omega(2);
    wz = omega(3);
    
    % assembles sensitivity matrix
    C = [h(wx,c)   0         0;
         0         h(wy,c)   0;
         0         0         h(wz,c)];
    
end



%==========================================================================
% h  Derivative of the saturation function.
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (1×1 double) independent variable
%   c       - (1×1 double) threshold value for saturation
%
% -------
% OUTPUT:
% -------
%   hx   	- (1×1 double) evaluation of h(x)
%
%==========================================================================
function hx = h(x,c)
    if iabs(x) < c
        hx = 1;
    else
        hx = 0;
    end
end