function func_err_time_indep(L_max,R_max,K_d)
% Function to check for errors in those functions that are
% time-independent (i.e. only dependent on K_d and not the kinetic rates)

% Throw an error message if L_max is not numeric or negative
if(~isnumeric(L_max) || L_max<=0)
   error('L_max should be a positive integer.');
end

% Throw an error message if R_max is not numeric or negative
if(~isnumeric(R_max) || R_max<=0)
   error('R_max should be a positive integer.');
end

% Throw an error message if K_d is not numeric or negative
if(~isnumeric(K_d) || K_d<=0)
   error('K_d should be a positive real number.');
end

end