% Script that implements a stochastic simulation of the 
% reversible heterodimerisation reaction for a comparative plot of 
% different binding and unbinding rates.

% Close any figures that are still open
close all

% Clear the workspace
clear

% Set the initial ligand copy-number
L_start = 10^3;
% Set the initial receptor copy-number
R_start = 10^1;
% Set the initial bound complex copy-number
B_start = 0;
% Set the total ligand and receptor copy-number
L_max = L_start + B_start;
R_max = R_start + B_start;
% Set the maximum bound complex copy-number
B_max = min(L_max,R_max);
% Set the maximum of the ligand and receptor copy-number
U_max = max(L_max,R_max);
% Set the minimum bound complex copy-number
B_min = 0;
% Calculate the number of states
n_states = length([L_start,R_start,B_start]);
% Set the start time of the simulation
t_start = 0;
% Set the end time of the simulation
t_end = 10^1;

% Get the standard matlab plotting colours
plot_col = get(gca,'ColorOrder');

% Set the text to be used on the y-axis of the panel rows
x1_ax_txt = {'low mean'; 'med entropy'; 'low ent rate'};
x2_ax_txt = {'med mean'; 'high entropy'; 'high ent rate'};
x3_ax_txt = {'high mean'; 'med entropy'; 'v high ent rate'};
x4_ax_txt = {'v high mean'; 'low entropy'; 'med ent rate'};

% Set the positions of the panels
sub_plot_pos = [0.125,0.725,0.45,0.225; 0.125,0.5,0.45,0.225;...
    0.125,0.275,0.45,0.225; 0.125,0.05,0.45,0.225;...
    0.575,0.725,0.45,0.225; 0.575,0.5,0.45,0.225;...
    0.575,0.275,0.45,0.225; 0.575,0.05,0.45,0.225];

% Set the unbinding rates of the columns
k_off_all = [1,10];
% Calculate the total number of columns
length_k_off_all = length(k_off_all);
% Set the binding rates of the rows
k_on_div_nu_all = [10.^(-4:-1);10.^(-3:0)]';
% Calculate the total number of rows
length_k_on_div_nu_all = length(k_on_div_nu_all);
% Set the letters of the panels
pan_let_init = 97;
pan_let = char(pan_let_init:pan_let_init+...
    (length_k_off_all*length_k_on_div_nu_all));

%%% Gillespie stochastic simulation algorithm %%%

% Loop through the unbinding rates
for j = 1:length_k_off_all;
    % Set the unbinding rate of the simulation
    k_off = k_off_all(j);
    % Loop through the binding rates
    for i = 1:length_k_on_div_nu_all
        % Set the binding rate of the simulation
        k_on_div_nu = k_on_div_nu_all(i,j);
        
        %%% Calculate analytical properties
        
        % Calculate the 2D dissociation constant
        K_d = k_off/k_on_div_nu;
        % Calculate the mean of the bound complex stationary distribution
        mean_B = func_mean(L_max,R_max,K_d);
        % Calculate the variance of the bound complex stationary distribution
        var_B = func_var(L_max,R_max,K_d);
        % Calculate the Shannon entropy of the bound complex stationary distribution
        ent_B = func_ent(L_max,R_max,K_d);
        % Calculate the mean reaction rate 
        % (i.e. the mean total number of reactions per second)
        reac_rate_B = 2*k_off*mean_B;
        % Calculate the entropy rate of the bound complex stationary distribution
        ent_rate_B = reac_rate_B*ent_B;
        % Calculate the EXPECTED total number of reactions
        reac_tot_exp_B = (t_end-t_start)*reac_rate_B;
        
        %%% Set initial values
        
        % Set the initial time until the next reaction
        tau_sim = t_start;
        % Set the initial time of the simulation
        t_sim = t_start;
        % Set the initial simulation step
        i_sim = 1;
        % Set the initial bound complex copy-number
        B_sim = B_start;
        % Set the state-change vector 
        % (increase bound complex copy-number for a binding reaction and 
        % decrease bound complex copy-number for an unbinding reaction)
        beta_sim = [1;-1];
        % Set a large storage matrix for the simulation output
        output_sim = zeros(10^6,2);
        
        % Start stochastic simulation
        while t_sim <= t_end
            % Calculate the propensity of a binding reaction
            alpha_sim(1) = k_on_div_nu*(B_max-B_sim)*(U_max-B_sim);
            % Calculate the propensity of an unbinding reaction
            alpha_sim(2) = k_off*B_sim;
            % Calculate the total reaction rate
            alpha_0_sim = sum(alpha_sim); 
            % Generate the time until next reaction from the exponential distribution
            tau_sim = exprnd(1/alpha_0_sim);
            % Record the duration of time spent in current state
            output_sim(i_sim,:) = [tau_sim,B_sim];
            % Generate the next reaction from the multinomial distribution
            reac_sim = mnrnd(1,alpha_sim/alpha_0_sim);
            % Update the system depending on which reaction occurred
            B_sim = B_sim + reac_sim*beta_sim;
            % Update the time of the simulation
            t_sim = t_sim + tau_sim;
            % Update the simulation step
            i_sim = i_sim + 1;
        end
        
        % Remove unused matrix storage
        output_sim(i_sim:length(output_sim),:) = [];
        
        % Calculate the SIMULATED total number of reactions
        reac_tot_sim_B = i_sim;
        
        %%% Print output to screen
        
        % Set text of the unbinding and binding rates of simulation
        disp_gen = [' (k_off = ', num2str(k_off),...
            ', k_on_div_nu = ', num2str(k_on_div_nu), ')'];
        % Set text of the expected total number of reactions
        disp_exp = ['Expected total number of reactions = ', ...
            num2str(round(reac_tot_exp_B,1))];
        % Ser text of the simulated total number of reactions
        disp_sim = ['Simulated total number of reactions = ', ...
            num2str(round(reac_tot_sim_B,1))];
        % Print the expected total number of reactions
        disp(strcat(disp_exp, disp_gen))
        % Print the simulated total number of reactions
        disp(strcat(disp_sim, disp_gen))
        
        %%% NOTE THAT THE DIFFERENCES BETWEEN THE EXPECTED AND SIMULATED 
        %%% TOTAL NUMBER OF REACTIONS ARE QUITE LARGE BECAUSE THE SIMULATION 
        %%% IS ONLY 10 SECONDS AS A DEFAULT. THESE ERRORS WILL REDUCE IF 
        %%% THE SIMULATION IS RUN FOR A LONGER PERIOD.
    
        %%% Plot output
        
        % Set panel
        subplot(length_k_on_div_nu_all,length_k_off_all,(2*i)-1 + (j-1));
        % Plot the time series of the bound complex copy-number
        plot(cumsum(output_sim(:,1)), output_sim(:,2),...
            'LineStyle', '-', 'Color', plot_col(1,:))
        % Plot a horizontal dashed line of the mean bound complex copy-number
        line([t_start,t_end],[mean_B,mean_B],...
            'LineStyle','--', 'Color', plot_col(2,:))
        % Plot a horizontal dotted line of the mean plus standard deviation
        line([t_start,t_end],[mean_B+sqrt(var_B),mean_B+sqrt(var_B)],...
            'LineStyle',':', 'Color', plot_col(2,:))
        % Plot a horizontal dotted line of the mean minus standard deviation
        line([t_start,t_end],[mean_B-sqrt(var_B),mean_B-sqrt(var_B)],...
            'LineStyle',':', 'Color', plot_col(2,:))
        % Set the limits of the x axes
        xlim([t_start,t_end])
        % Set the limits of the y axes
        ylim([B_min, B_max])
        % Set the panel letter
        text(-0.12*(t_end-t_start),0.95*(B_max-B_min),...
            strcat('\fontsize{11} \bf', {' '}, pan_let((2*i)-1 + (j-1))))
        % Get current axes handle
        ax = gca;
        % Remove tick marks on x axes
        ax.XTick = [];
        % Remove tick marks on y axes
        ax.YTick = [];
        % Set the outer position of the axes
        ax.OuterPosition = sub_plot_pos(i + (length_k_on_div_nu_all*(j-1)),:);
        % Add text of the binding rate to each panel and adapt the box
        % surrounding the binding rate depending on whether it overlaps
        % with the stochastic simulation.
        if(i==1 || i==2)
            te = text(0.45,8,strcat('$k_\textrm{on}/\nu$=',...
                num2str(k_on_div_nu_all(i,j)),'/s'),'Interpreter', 'latex');
            if(i==2)
                te.EdgeColor = 'k';
                te.BackgroundColor = 'w';
                te.Margin = 1;
            end
        else
            te = text(0.45,2,strcat('$k_\textrm{on}/\nu$=',...
                num2str(k_on_div_nu_all(i,j)),'/s'),'Interpreter', 'latex');
        end
        % Add text of the unbinding rate and the text to the left of each 
        % row of panels
        if(i==1 && j==1)
            text(4,12,strcat('$k_\textrm{off}$=',...
                num2str(k_off_all(j)),'/s'), 'Interpreter', 'latex')
            text(-4.25,5,x1_ax_txt, 'Interpreter', 'latex')
        end
        if(i==2 && j==1)
            text(-4.25,5,x2_ax_txt, 'Interpreter', 'latex')
        end
        if(i==3 && j==1)
            text(-4.25,5,x3_ax_txt, 'Interpreter', 'latex')
        end
        if(i==4 && j==1)
            text(-4.25,5,x4_ax_txt, 'Interpreter', 'latex')
        end
        if(i==1 && j==2)
            text(3.5,12,strcat('$k_\textrm{off}$=',...
                num2str(k_off_all(j)),'/s'), 'Interpreter', 'latex')
        end
        hold on
    end
end
hold off
% Add text indicating the duration of the simulation
text(-3,-2,strcat('Time (', num2str(t_start), '$\rightarrow$',...
    num2str(t_end), ' s)'), 'Interpreter', 'latex')

%%% Write figure to file

% Save as a .eps file 
print('fig_stoc_sim', '-depsc')
% Save as a .emf file
print('fig_stoc_sim', '-dmeta')
