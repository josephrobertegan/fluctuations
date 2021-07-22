function func_err_time_dep(L_max,R_max,K_d,k_off)
% Function to check for errors in those functions that are
% time-dependent (i.e. dependent on binding rates)

% Throw an error message if L_max is not numeric or not positive
if(~isnumeric(L_max) || L_max<=0)
   error('L_max should be a positive integer.');
end

% Throw an error message if R_max is not numeric or not positive
if(~isnumeric(R_max) || R_max<=0)
   error('R_max should be a positive integer.');
end

% Throw an error message if K_d is not numeric or not positive
if(~isnumeric(K_d) || K_d<=0)
   error('K_d should be a positive real number.');
end

% Throw an error message if k_off is not numeric or not positive
if(~isnumeric(k_off) || k_off<=0)
   error('k_off should be a positive real number.');
end

end