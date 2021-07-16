function result = func_ent(L_max,R_max,K_d)
% Computational calculation for the Shannon entropy of the bound complex copy-number 
% stationary distribution of the reversible heterodimerisation reaction.

% Calculate stationary distribution
stat_dist = func_stat_dist(L_max,R_max,K_d);
% Calculate the Shannon entropy
result = -log2(prod(stat_dist.^stat_dist));
end