function result = func_ent_rate(L_max,R_max,K_d,k_off)
% Computational calculation for the entropy rate of the bound complex copy-number 
% stationary distribution of the reversible heterodimerisation reaction.

% Check that user inputs are positive
func_err_time_dep(L_max,R_max,K_d,k_off)

% Calculate the mean
mean = func_mean(L_max,R_max,K_d);
% Calculate the Shannon entropy
entropy = func_ent(L_max,R_max,K_d);
% Calculate the entropy rate
result = 2*k_off*mean*entropy;
end