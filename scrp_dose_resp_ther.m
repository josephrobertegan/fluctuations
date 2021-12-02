% Script that produces a figure of the theoretical dose-response plots 
% where response is the variance rate, entropy rate, or mean signalling 
% rate and dose is the ligand copy-number. The variance rate and entropy 
% rate are calculated analytically, and the mean signalling rate is 
% calculated via stochastic simulations.

% Close any figures that are still open
close all

% Clear the workspace
clear

% Set constants for positioning panel letters
pan_let_x_scale = -0.1;
pan_let_y_scale = 1.15;

% Set constants for positioning the sub-plots
plot_sub_x_left = 0.09;
plot_sub_x_right = 0.59;
plot_sub_y_bottom = 0.1;
plot_sub_y_top = 0.59;

% Set constants for the size of the sub-plots
plot_sub_width = 0.4;
plot_sub_height = 0.345;

%%% General parameters required for calculating dose response

% Set all of the unique receptor copy-numbers
R_max_all = [10^1,10^2];
% Calculate the total number of unique receptor copy-numbers
R_max_all_length = length(R_max_all);
% Set the minimum ligand copy-number on a log10 scale
L_max_lower = 0;
% Set the maximum ligand copy-number on a log10 scale
L_max_upper = 3;
% Set the increment for the ligand copy-number on a log10 scale
L_max_inc = 0.01;
% Calculate all of the unique ligand copy-numbers
L_max_all = unique(round(10.^(L_max_lower:L_max_inc:L_max_upper)));
% Calculate the total number of unique ligand copy-numbers
L_max_all_length = length(L_max_all);
% Calculate all of the unique ligand copy-numbers on a log10 scale
log10_L_max_all = log10(L_max_all);

% Set the fixed unbinding rate
k_off = 10^0;
% Set the minimum binding rate on a log10 scale
k_on_div_nu_lower = -4;
% Set the maximum binding rate on a log10 scale
k_on_div_nu_upper = -1;
% Set the increment for the binding rate on a log10 scale
k_on_div_nu_inc = 1;
% Calculate all of the binding rates
k_on_div_nu_all = 10.^(-(-k_on_div_nu_upper:k_on_div_nu_inc:-k_on_div_nu_lower));
% Calculate the total number of binding rates
k_on_div_nu_all_length = length(k_on_div_nu_all);
% Calculate all of the 2D dissociation constants
K_d_all = k_off./k_on_div_nu_all;
% Calculate the total number of 2D dissociation constants
K_d_all_length = length(K_d_all);

% Set the legend text for the plots
leg_txt = strsplit(num2str(log10(K_d_all)));
% Set the legend ttitle for the plots
leg_tit = '$\log_{10}(K_\textrm{d})$';

% Set a scaling for extending the axes by a small fraction
ax_scale = 0.025;
% Calculate the minimum ligand copy-number on a log10 scale
log10_x_min = log10(min(L_max_all));
% Calculate the maximum ligand copy-number on a log10 scale
log10_x_max = log10(max(L_max_all));
% Calculate the extension at both ends of the x axis
x_ext = ax_scale*(log10_x_max-log10_x_min);
% Calculate a lower bound for the x axis
x_lower = log10_x_min - x_ext;
% Calculate an upper bound for the x axis
x_upper =  log10_x_max + x_ext;

%%% Parameters required for calculating signalling rate via stochastic
%%% simulation

% Set the initial inactive receptor copy-number
R_Istart = R_max_all(1);
% Set the initial active receptor copy-number
R_Astart = 0;
% Set the initial bound complex copy-number
B_start = 0;
% Set the initial signal copy-number
S_start = 0;

% Set the number of repetitions of the simulation
rep_length = 10;
% Set the start time of the simulation
t_start = 0;
% Set the time of the simulation at which to terminate the simulation
t_end = 10^4;
% Set the size of the signal copy-number at which to terminate the simulation
S_end = 10^4;

% Set a storage array for the variance rate of the bound complex stationary distribution
var_rate_B = zeros(L_max_all_length, k_on_div_nu_all_length, R_max_all_length);
% Set a storage array for the entropy rate of the bound complex stationary distribution
ent_rate_B = zeros(L_max_all_length, k_on_div_nu_all_length, R_max_all_length);
% Set a storage array for the signalling rate
sig_rate_B = zeros(L_max_all_length,k_on_div_nu_all_length,rep_length);
% Set a storage array for the mean signalling rate
mean_sig_rate_tmp = zeros(L_max_all_length,k_on_div_nu_all_length);

% Loop through all values of the ligand copy-number
for i = 1:L_max_all_length
    % Set the ligand copy-number of current loop
    L_max = L_max_all(i);
    % Loop through all values of the unbinding rate
    for j = 1:k_on_div_nu_all_length
        % Set the unbinding rate of current loop
        k_on_div_nu = k_on_div_nu_all(j);
        % Calculate the 2D dissociation constant
        K_d = k_off/k_on_div_nu;
        % Loop through all values of the receptor copy-number
        for k = 1:R_max_all_length
            % Set the receptor copy-number of current loop
            R_max = R_max_all(k);
            % Calculate the variance rate of the bound complex stationary distribution
            var_rate_B(i,j,k) = func_var_rate(L_max,R_max,K_d,k_off);
            % Calculate the entropy rate of the bound complex stationary distribution
            ent_rate_B(i,j,k) = func_ent_rate(L_max,R_max,K_d,k_off);
        end
    end
end

% Loop through all values of the stochastic simulation repetitions
for k=1:rep_length
    % Set the stochastic simulation repetition of current loop
    rep = k;
    % Loop through all values of the ligand copy-number
    for i=1:L_max_all_length
        % Set the ligand copy-number of current loop
        L_start = L_max_all(i);
        % Loop through all values of the unbinding rate
        for j=1:k_on_div_nu_all_length
            % Set the unbinding rate of current loop
            k_on_div_nu = k_on_div_nu_all(j);
            % Set the initial time until the next reaction
            tau_sim = t_start;
            % Set the initial time of the simulation
            t_sim = t_start;
            % Set the initial simulation step 
            i_sim = 1;
            % Set the initial ligand copy-number
            L_sim = L_start;
            % Set the initial inactive receptor copy-number
            R_Isim = R_Istart;
            % Set the initial active receptor copy-number
            R_Asim = R_Astart;
            % Set the initial bound complex copy-number
            B_sim = B_start;
            % Set the initial signal copy-number
            S_sim = S_start;
            % Set the initial state vector
            X_sim = [L_sim, R_Isim, R_Asim, B_sim, S_sim];
            % Set the state-change vector corresponding to the propensities 
            % (i.e. corresponding to the alpha_sim values below)
            beta_sim = [-1,-1,0,1,0; -1,0,-1,1,0; 1,0,1,-1,0; 0,1,-1,0,0; 0,0,0,0,1];
            % Calculate the total number of states
            n_states = length(X_sim);
            % Set a large storage matrix for the simulation output
            output_sim = zeros(10^6,n_states+1);
            % Start the stochastic simulation and continue running it 
            % until either the signal copy-number or the time of the 
            % simulation has reached its limit.
            while (X_sim(n_states) <= S_end) && (t_sim <= t_end)
                % Calculate the propensity of binding of L and R_I to form B
                alpha_sim(1) = k_on_div_nu*X_sim(1)*X_sim(2);
                % Calculate the propensity of binding of L and R_A to form B
                alpha_sim(2) = k_on_div_nu*X_sim(1)*X_sim(3);
                % Calculate the propensity of unbinding of B to form L and R_A
                alpha_sim(3) = k_off*X_sim(4);
                % Calculate the propensity of reverting of R_A to R_I
                alpha_sim(4) = k_off*X_sim(3);
                % Calculate the propensity of signal generation of S
                alpha_sim(5) = k_off*X_sim(3)*X_sim(4);
                % Calculate the total reaction rate
                alpha_0_sim = sum(alpha_sim); 
                % Generate the time until next reaction from exponential distribution
                tau_sim = exprnd(1/alpha_0_sim);
                % Record the duration of time spent in the current state
                output_sim(i_sim,:) = [tau_sim,X_sim];
                % Generate the next reaction from multinomial distribution
                reac_sim = mnrnd(1,alpha_sim/alpha_0_sim);
                % Update the system depending on which event occurred
                X_sim = X_sim + reac_sim*beta_sim;
                % Update the time of the simulation
                t_sim = t_sim + tau_sim;
                % Update the simulation step
                i_sim = i_sim + 1;
            end
            
            % Remove unused matrix storage
            output_sim(i_sim:length(output_sim),:) = [];
            
            % Calculate the signalling rate
            sig_rate_B(i,j,k) = X_sim(n_states)/t_sim;
            
            % Set text of the repetition number        
            disp_rep = strcat('Rep=', num2str(rep), ',');
            % Set text of the ligand copy-number on a log10 scale
            disp_L_max = strcat('Lmax=', num2str(round(log10(L_start),...
                -log10(L_max_inc))),',');
            % Set text of the 2D dissociation constant on a log10 scale
            disp_K_d = strcat('Kd=', num2str(round(log10(k_off/k_on_div_nu),...
                -log10(k_on_div_nu_inc))), ',');
            % Set text of the signalling rate
            disp_S = strcat('SigRate=', num2str(round(sig_rate_B(i,j,k),2)));
            % Print text above to the screen to keep track of simulations
            disp(strcat(disp_rep,disp_L_max,disp_K_d,disp_S))
        end
    end
end

% Calculate mean signalling rate by taking the average of the signalling
% rate across all of the repetitions of the stochastic simulations
for k = 1:rep_length
    mean_sig_rate_tmp = mean_sig_rate_tmp + sig_rate_B(:,:,k);
end
mean_sig_rate_B = mean_sig_rate_tmp/rep_length;

% Set the first subplot
subplot(2,2,1)
% Plot the entropy rate with smaller R_max against (log10) ligand copy-number 
plot(log10_L_max_all,ent_rate_B(:,:,1))
% Calculate the minimum entropy rate over all unbinding rates and ligand
% copy-numbers
y_min_1 = min(min(ent_rate_B(:,:,1)));
% Calculate the maximum entropy rate over all unbinding rates and ligand
% copy-numbers
y_max_1 = max(max(ent_rate_B(:,:,1)));
% Calculate the extension to add to the y-axis
y_ext_1 = ax_scale*(y_max_1-y_min_1);
% Calculate the lower y axis value
y_lower_1 = y_min_1-y_ext_1;
% Calculate the upper y axis value
y_upper_1 =  y_max_1+y_ext_1;
% Set the y axis limits
ylim([y_lower_1,y_upper_1])
% Set the x axis limits
xlim([x_lower,x_upper])
% Add the panel letter
text(x_lower + pan_let_x_scale*(x_upper-x_lower),...
    y_lower_1 + pan_let_y_scale*(y_upper_1-y_lower_1),'\fontsize{12} \bf a')
% Get the current axis handle
ax = gca;
% Set the position of the axes
ax.Position = [plot_sub_x_left, plot_sub_y_top, plot_sub_width, plot_sub_height];
% Add the legend text
l = legend(leg_txt, 'Location', 'northwest', 'Interpreter', 'latex');
% Add the legend title
title(l, leg_tit, 'Interpreter', 'latex')
% Remove the box surrounding the legend
l.Box = 'off';
% Add the plot title
title(strcat('$R_{\textrm{max}}=',num2str(R_max_all(1)),'$'), 'Interpreter', 'latex')
% Add the x-axis label
xlabel('$\log_{10}(L_{\textrm{max}})$', 'Interpreter', 'latex')
% Add the y-axis label
ylabel('Entropy rate (bits/s)', 'Interpreter', 'latex')

% Set the second subplot
subplot(2,2,2)
% Plot variance rate with smaller R_max against (log10) ligand copy-number 
plot(log10_L_max_all,var_rate_B(:,:,1))
% Calculate minimum variance rate over all unbinding rates and ligand
% copy-numbers
y_min_2 = min(min(var_rate_B(:,:,1)));
% Calculate maximum variance rate over all unbinding rates and ligand
% copy-numbers
y_max_2 = max(max(var_rate_B(:,:,1)));
% Calculate the extension to add to the y-axis
y_ext_2 = ax_scale*(y_max_2-y_min_2);
% Calculate the lower y axis value
y_lower_2 = y_min_2-y_ext_2;
% Calculate the upper y axis value
y_upper_2 =  y_max_2+y_ext_2;
% Set the y axis limits
ylim([y_lower_2,y_upper_2])
% Set the x axis limits
xlim([x_lower,x_upper])
% Add the panel letter
text(x_lower + pan_let_x_scale*(x_upper-x_lower),...
    y_lower_2 + pan_let_y_scale*(y_upper_2-y_lower_2),'\fontsize{12} \bf b')
% Get the current axis handle
ax = gca;
% Set the position of the axes
ax.Position = [plot_sub_x_right, plot_sub_y_top, plot_sub_width, plot_sub_height];
% Add the legend text
l = legend(leg_txt, 'Location', 'northwest', 'Interpreter', 'latex');
% Add the legend title
title(l, leg_tit, 'Interpreter', 'latex')
% Remove the box surrounding the legend
l.Box = 'off';
% Add the plot title
title(strcat('$R_{\textrm{max}}=',num2str(R_max_all(1)),'$'), 'Interpreter', 'latex')
% Add the x-axis label
xlabel('$\log_{10}(L_{\textrm{max}})$', 'Interpreter', 'latex')
% Add the y-axis label
ylabel('Variance rate (s$^{-1}$)', 'Interpreter', 'latex')

% Set the third subplot
subplot(2,2,3)
% Plot entropy rate with larger R_max against (log10) ligand copy-number
plot(log10_L_max_all,ent_rate_B(:,:,2))
% Calculate minimum entropy rate over all unbinding rates and ligand
% copy-numbers
y_min_3 = min(min(ent_rate_B(:,:,2)));
% Calculate maximum entropy rate over all unbinding rates and ligand
% copy-numbers
y_max_3 = max(max(ent_rate_B(:,:,2)));
% Calculate the extension to add to the y-axis
y_ext_3 = ax_scale*(y_max_3-y_min_3);
% Calculate the lower y axis value
y_lower_3 = y_min_3-y_ext_3;
% Calculate the upper y axis value
y_upper_3 =  y_max_3+y_ext_3;
% Set the y axis limits
ylim([y_lower_3,y_upper_3])
% Set the x axis limits
xlim([x_lower,x_upper])
% Add the panel letter
text(x_lower + pan_let_x_scale*(x_upper-x_lower),...
    y_lower_3 + pan_let_y_scale*(y_upper_3-y_lower_3),'\fontsize{12} \bf c')
% Get the current axis handle
ax = gca;
% Set the position of the axes
ax.Position = [plot_sub_x_left, plot_sub_y_bottom, plot_sub_width, plot_sub_height];
% Add the legend text
l = legend(leg_txt, 'Location', 'northwest', 'Interpreter', 'latex');
% Add the legend title
title(l, leg_tit, 'Interpreter', 'latex')
% Remove the box surrounding the legend
l.Box = 'off';
% Add the plot title
title(strcat('$R_{\textrm{max}}=',num2str(R_max_all(2)),'$'), 'Interpreter', 'latex')
% Add the x-axis label
xlabel('$\log_{10}(L_{\textrm{max}})$', 'Interpreter', 'latex')
% Add the y-axis label
ylabel('Entropy rate (bits/s)', 'Interpreter', 'latex')

% Set the fourth subplot
subplot(2,2,4)
% Plot mean signalling rate with smaller R_max against (log10) ligand copy-number
plot(log10_L_max_all, mean_sig_rate_B)
% Calculate minimum mean signalling rate over all unbinding rates and ligand
% copy-numbers
y_min_4 = min(min(mean_sig_rate_B));
% Calculate maximum mean signalling rate over all unbinding rates and ligand
% copy-numbers
y_max_4 = max(max(mean_sig_rate_B));
% Calculate the extension to add to the y-axis
y_ext_4 = ax_scale*(y_max_4 - y_min_4);
% Calculate the lower y axis value
y_lower_4 = y_min_4 - y_ext_4;
% Calculate the upper y axis value
y_upper_4 =  y_max_4 + y_ext_4;
% Set the y axis limits
ylim([y_lower_4,y_upper_4])
% Set the x axis limits
xlim([x_lower,x_upper])
% Add the panel letter
text(x_lower + pan_let_x_scale*(x_upper-x_lower),...
    y_lower_4 + pan_let_y_scale*(y_upper_4-y_lower_4),'\fontsize{12} \bf d')
% Get the current axis handle
ax = gca;
% Set the position of the axes
ax.Position = [plot_sub_x_right, plot_sub_y_bottom, plot_sub_width, plot_sub_height];
% Add the legend text
l = legend(leg_txt, 'Location', 'northwest', 'Interpreter', 'latex');
% Add the legend title
title(l, leg_tit, 'Interpreter', 'latex')
% Remove the box surrounding the legend
l.Box = 'off';
% Add the plot title
title(strcat('$R_{\textrm{max}}=',num2str(R_Istart),'$'), 'Interpreter', 'latex')
% Add the x-axis label
xlabel('$\log_{10}(L_{\textrm{max}})$', 'Interpreter', 'latex')
% Add the y-axis label
ylabel('Signaling rate (s$^{-1}$)', 'Interpreter', 'latex')

%%% Write figure to file

% Save as a .eps file
print('fig_dose_resp_ther', '-depsc')
% Save as a .emf file
print('fig_dose_resp_ther', '-dmeta')
