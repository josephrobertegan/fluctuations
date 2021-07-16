% Script that gives examples of how to use the various functions 

% Clear the workspace
clear

% Example values of L_max, R_max, k_off and k_on_div_nu
Lmax = 10; % Dimensionless
Rmax = 10; % Dimensionless
koff = 1; % Dimension of time^{-1}
kon_div_nu = 0.1; % Dimension of time^{-1}

% Calculate example 2D K_d
Kd = koff/kon_div_nu; % Dimensionless

% Calculate example stationary distribution of the bound complex
% copy-number
stat_exam = func_stat_dist(Lmax,Rmax,Kd);

% Calculate example mean of stationary distribution
mean_exam = func_mean(Lmax,Rmax,Kd);

% Calculate example variance of stationary distribution
var_exam = func_var(Lmax,Rmax,Kd);

% Calculate example entropy of stationary distribution
ent_exam = func_ent(Lmax,Rmax,Kd);

% Calculate example variance rate of stationary distribution
var_rate_exam = func_var_rate(Lmax,Rmax,Kd,koff);

% Calculate example entropy rate of stationary distribution
ent_rate_exam = func_ent_rate(Lmax,Rmax,Kd,koff);

% Create a temporary function
func_var_rate_opt_tmp = @(K_d) func_var_rate_opt(Lmax,Rmax,K_d);
% Calculate the value of K_d that maximises the variance rate for fixed 
% example values of L_max and R_max
K_d_opt_var_rate = fzero(func_var_rate_opt_tmp, max(Lmax,Rmax)/2);

% Create a temporary function
func_ent_rate_opt_tmp = @(K_d) func_ent_rate_opt(Lmax,Rmax,K_d);
% Calculate the value of K_d that maximises the entropy rate for fixed 
% example values of L_max and R_max
K_d_opt_ent_rate = fzero(func_ent_rate_opt_tmp, max(Lmax,Rmax)/4);
