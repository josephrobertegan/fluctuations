function result = func_var(L_max,R_max,K_d)
% Computational calculation for the variance of the bound complex copy-number 
% stationary distribution of the reversible heterodimerisation reaction.

% Calculate the mean
mean = func_mean(L_max,R_max,K_d);
% Calculate the variance
result = mean*(L_max + R_max + K_d - mean) - L_max*R_max;
end