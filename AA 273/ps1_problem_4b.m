%% ps1_problem_4b
% Problem Set 1, Problem 4b
% AA 273 - State Estimation and Filtering for Robotic Perception
%
% Author: Tamas Kis
% Last Update: 2021-08-18



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;



%% SOLUTION

% initial condition
pi_0 = [1;0];

% transition matrix
T = [0.999,0.05;0.001,0.95];

% ones vector
ones_vector = [1;1];

% pi_(t-1) for t=1 (i.e. initial condition)
pi_tm1 = pi_0;

% histogram filter
for t = 1:99
    
    % determines yt based on day
    if t ~= 99
        yt = 0;
    else
        yt = 1;
    end
    
    % measurement likelihood matrix
    if yt == 0
        M = [0.9,0;0,0.01];
    else
        M = [0.1,0;0,0.99];
    end
    
    % new belief vector
    pi_t = (M*T*pi_tm1)/((M*ones_vector)'*T*pi_tm1);
    
    % storers belief vector for next iteration
    pi_tm1 = pi_t;
    
end

% probability they have illness given test history (second element of 
% belief vector)
pi_t(2)