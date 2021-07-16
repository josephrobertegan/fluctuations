function result = func_ent_rate_opt(L_max,R_max,K_d)
% Function to allow for the numerical calculation of the optimal K_d that 
% maximises the entropy rate for fixed L_max and R_max.

% Calculate the mean
mean = func_mean(L_max,R_max,K_d);
% Calculate the variance
variance = func_var(L_max,R_max,K_d);
% Calculate the Shannon entropy
entropy = func_ent(L_max,R_max,K_d);
% Calculate all possible values that the bound complex copy-number can take
B_all = 0:min(L_max,R_max);
% Calculate stationary distribution
stat_dist = func_stat_dist(L_max,R_max,K_d);
% Calculate J as given in Supplementary Information of paper
J = -log2(prod(stat_dist.^(B_all.*stat_dist)));
% Equate the following result to zero and solve for K_d
result = (mean^2 - variance)*entropy - mean*J;
end