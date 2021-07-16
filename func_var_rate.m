function result = func_var_rate(L_max,R_max,K_d,k_off)
% Computational calculation for the variance rate of the bound complex copy-number 
% stationary distribution of the reversible heterodimerisation reaction.

% Calculate the mean
mean = func_mean(L_max,R_max,K_d);
% Calculate the variance
variance = func_var(L_max,R_max,K_d);
% Calculate the variance rate
result = 2*k_off*mean*variance;
end