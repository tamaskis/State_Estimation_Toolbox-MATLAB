%-------------------------------------------------------------------------%

% TAMAS KIS

% AA 273 - State Estimation and Filtering for Robotic Perception
% Problem Set 1 - Problem 3a

%-------------------------------------------------------------------------%



%% SCRIPT SETUP

% clears variables and command window, closes all figures
clear;
clc;
close all;

% adds path to Estimation Toolbox
addpath("..");



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