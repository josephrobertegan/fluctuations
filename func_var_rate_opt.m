function result = func_var_rate_opt(L_max,R_max,K_d)
% Function to allow for the numerical calculation of the optimal K_d that 
% maximises the variance rate for fixed L_max and R_max.

% Calculate the mean
mean = func_mean(L_max,R_max,K_d);
% Calculate the variance
variance = func_var(L_max,R_max,K_d);
% Equate the following result to zero and solve for K_d
result = K_d*(mean^2) - variance*(variance + mean*(L_max + R_max + K_d - 2*mean));
end