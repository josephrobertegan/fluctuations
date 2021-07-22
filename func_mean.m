function result = func_mean(L_max,R_max,K_d)
% Computational calculation for the mean of the bound complex copy-number 
% stationary distribution of the reversible heterodimerisation reaction.

% Check that user inputs are positive
func_err_time_indep(L_max,R_max,K_d)

% Calculate all possible values that the bound complex copy-number can take
B_all = 0:min(L_max,R_max);
% Calculate stationary distribution
stat_dist = func_stat_dist(L_max,R_max,K_d);
% Calculate expectation of stationary distribution (i.e. the mean)
result = B_all*stat_dist';
end
