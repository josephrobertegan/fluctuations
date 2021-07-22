function result = func_var_rate(L_max,R_max,K_d,k_off)
% Computational calculation for the variance rate of the bound complex copy-number 
% stationary distribution of the reversible heterodimerisation reaction.

% Check that user inputs are positive
func_err_time_dep(L_max,R_max,K_d,k_off)

% Calculate the mean
mean = func_mean(L_max,R_max,K_d);
% Calculate the variance
variance = func_var(L_max,R_max,K_d);
% Calculate the variance rate
result = 2*k_off*mean*variance;
end