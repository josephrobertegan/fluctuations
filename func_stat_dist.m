function result = func_stat_dist(L_max,R_max,K_d)
% Computational calculation for the stochastic stationary solution of  the 
% bound complex copy-number of the reversible heterodimerisation reaction.
% This approach avoids computational overflow issues that can occur with 
% large values of B_max (i.e. it avoids indeterminate forms of the type
% inf/inf)

% Calculate B_max and U_max
B_max = min(L_max,R_max);
U_max = max(L_max,R_max);
% Calculate all possible values that the bound complex copy-number can take
B_all = 0:B_max;
% Calculate the negative of the natural log of a(B) where p(B) = a(B)/sum(a(B))
neg_log_a_B = B_all*log(K_d) - gammaln(B_max+1) - gammaln(U_max+1) + ...
    gammaln(B_all+1) + gammaln(B_max-B_all+1) + gammaln(U_max-B_all+1);
% Calculate the minimum of the value above
min_neg_log_a_B = min(neg_log_a_B);
% Calculate a scaled value of a(B)
a_B_scaled = exp(min_neg_log_a_B-neg_log_a_B);
% Calculate p(B) in which the scaling of a(B) cancels in the ratio
result = a_B_scaled/sum(a_B_scaled);
end
