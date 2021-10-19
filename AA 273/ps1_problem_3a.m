%% ps1_problem_3a
% Problem Set 1, Problem 3a
% AA 273 - State Estimation and Filtering for Robotic Perception
%
% Author: Tamas Kis
% Last Update: 2021-08-18



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;



%% SOLUTION

% initial prior
prior = 0.99*0.001/0.10089;

% recursive Bayesian estimation
for n = 1:1
    posterior = (0.99*prior)/(0.1*(1-prior)+(0.99*prior));
    prior = posterior;
end

% final posterior
posterior