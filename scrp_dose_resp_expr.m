% Script that produces a figure of the experimental dose-response plots
% where response is some read-out of T cell activation and dose is the 
% ligand copy-number.

% Close any figures that are still open
close all

% Clear the workspace
clear

% Set the constants for positioning the panel letters
pan_let_x_scale = -0.1;
pan_let_y_scale = 1.15;

% Set the constants for positioning the sub-plots
plot_sub_x_left = 0.1;
plot_sub_x_right = 0.6;
plot_sub_y_bottom = 0.1;
plot_sub_y_top = 0.6;

% Set the constants for the size of the sub-plots
plot_sub_width = 0.395;
plot_sub_height = 0.325;

% Set a scaling for extending the axes by a small fraction
ax_scale = 0.05;

% Set a constant for converting molar concentrations to micro-molar concentrations
molar2micromolar = 10^6;
% Set a constant for converting pg/ml to ng/ml
pg2ng = 10^(-3);
% Set a constant for converting per minute to per second
permin2persec = 60^(-1);

% Set a label for the x axes
x_lab = '$\log_{10}([L_\textrm{max}])$';

% Set a title for the legends
leg_tit = '$\log_{10}([K_\textrm{d}])$';

% Get the standard matlab plotting colours
plot_col = get(gca,'ColorOrder');

%%% Plot andersen2001 data

% Set the 3D Kd values taken from andersen2001
K_d_3D_Andersen = molar2micromolar*[2.3*10^(-9), 2.8*10^(-8), 7.7*10^(-8),...
    1.3*10^(-7), 1.7*10^(-7), 3.2*10^(-7), 6.2*10^(-6)];
% Set the legend text
leg_txt_Andersen = {num2str(round(log10(K_d_3D_Andersen),1)')};

% Set the name of the file containing the data
file_name_Andersen = 'data_andersen_2001.csv';
% Read the data
data_Andersen = readtable(file_name_Andersen, 'ReadVariableNames', 1);

% Calculate the total number of different ligands
n_ligands_Andersen = length(data_Andersen.Properties.VariableNames) - 1;
% Calculate the dose values on a log10 scale
dose_Andersen = log10(data_Andersen.(char(data_Andersen.Properties.VariableNames(1))));
% Calculate the total number of dose values
n_dose_Andersen = length(dose_Andersen);

% Calculate the maximum dose
x_max_Andersen = max(dose_Andersen);
% Calculate the minimum dose
x_min_Andersen = min(dose_Andersen);
% Calculate the extension to add to the x-axis
x_ext_Andersen = ax_scale*(x_max_Andersen-x_min_Andersen);
% Calculate the lower x axis value
x_lower_Andersen = -1.5; % Hard-coded to fit legend in
% Calculate the upper x axis value
x_upper_Andersen = x_max_Andersen + x_ext_Andersen;

% Get all of the response data
response_data_Andersen = table2array(data_Andersen(:,2:end));
% Calculate the minimum response data over all 2D dissociation constants
% (i.e. different ligands) and ligand copy-numbers (i.e. different doses)
y_min_Andersen = min(min(response_data_Andersen));
% Calculate the maximum response data over all 2D dissociation constants
% (i.e. different ligands) and ligand copy-numbers (i.e. different doses)
y_max_Andersen = max(max(response_data_Andersen));
% Calculate the extension to add to the y-axis
y_ext_Andersen = ax_scale*(y_max_Andersen-y_min_Andersen);
% Calculate the lower y axis value
y_lower_Andersen = y_min_Andersen - y_ext_Andersen;
% Calculate the lower y axis value
y_upper_Andersen = y_max_Andersen + y_ext_Andersen;

% Set a storage matrix for analysing the response data
response_Andersen = zeros(n_dose_Andersen, n_ligands_Andersen);
% Get response data from first ligand
response_Andersen(:,1) = data_Andersen.(char(data_Andersen.Properties.VariableNames(2)));
% Calculate a vector of NaN values where 0 represents NaN and 1 otherwise
index_nan_Andersen = ~isnan(response_Andersen(:,1));
% Set the first subplot
subplot(2,2,1)
% Plot dose against response for first ligand (for all non NaN values)  
plot(dose_Andersen(index_nan_Andersen), response_Andersen(index_nan_Andersen,1),...
     'Marker', '*', 'Color', plot_col(1,:));
% Loop through remaining ligands and plot dose response as above
for i=1:(n_ligands_Andersen-1)
    response_Andersen(:,i+1) = data_Andersen.(char(data_Andersen.Properties.VariableNames(i+2)));
    index_nan_Andersen = ~isnan(response_Andersen(:,i+1));
    line(dose_Andersen(index_nan_Andersen), response_Andersen(index_nan_Andersen,i+1),...
        'Marker', '*', 'Color', plot_col(i+1,:))
end
% Set the x axis limits
xlim([x_lower_Andersen, x_upper_Andersen])
% Set the y axis limits
ylim([y_lower_Andersen,y_upper_Andersen])
% Add the panel letter
text(x_lower_Andersen + pan_let_x_scale*(x_upper_Andersen-x_lower_Andersen),...
    y_lower_Andersen + pan_let_y_scale*(y_upper_Andersen-y_lower_Andersen),...
    '\fontsize{12} \bf a')
% Get the current axis handle
ax = gca;
% Set the position of the axes
ax.Position = [plot_sub_x_left, plot_sub_y_top, plot_sub_width, plot_sub_height];
% Add the legend text
l = legend(leg_txt_Andersen, 'Interpreter', 'latex', 'Location', 'best');
% Add the legend title
title(l, leg_tit, 'Interpreter', 'latex')
% Remove the box surrounding the legend
legend('boxoff')
% Add the plot title
title('Andersen $\emph{et al.}$ 2001', 'Interpreter', 'latex')
% Add the x-axis label
xlabel(x_lab, 'Interpreter', 'latex')
% Add the y-axis label
ylabel('NFAT activated T cells (\%)', 'Interpreter', 'latex')

%%% Plot chmielewski2004 data

% Set 3D Kd values taken from Chmielewski2004
K_d_3D_Chmielewski = molar2micromolar*[1.5*10^(-11), 1.2*10^(-10), 1*10^(-9),...
    1.6*10^(-8), 3.2*10^(-7)];
% Set the legend text
leg_txt_Chmielewski = {num2str(round(log10(K_d_3D_Chmielewski),1)')};

% Set the name of the file containing the data
file_name_Chmielewski = 'data_chmielewski_2004.csv';
% Read the data
data_Chmielewski = readtable(file_name_Chmielewski, 'ReadVariableNames', 1);

% Calculate the total number of different ligands
n_ligands_Chmielewski = length(data_Chmielewski.Properties.VariableNames) - 1;
% Calculate the dose values of a log10 scale
dose_Chmielewski = log10(data_Chmielewski.(char(data_Chmielewski.Properties.VariableNames(1))));
% Calculate the total number of dose values
n_dose_Chmielewski = length(dose_Chmielewski);

% Calculate the maximum dose
x_max_Chmielewski = max(dose_Chmielewski);
% Calculate the minimum dose
x_min_Chmielewski = min(dose_Chmielewski);
% Calculate the extension to add to the x-axis
x_ext_Chmielewski = ax_scale*(x_max_Chmielewski-x_min_Chmielewski);
% Calculate the lower x-axis value
x_lower_Chmielewski = -2; % Hard-coded to fit legend in
% Calculate the upper x-axis value
x_upper_Chmielewski = x_max_Chmielewski + x_ext_Chmielewski;

% Get all of the response data
response_data_Chmielewski = table2array(data_Chmielewski(:,2:end));
% Calculate the minimum response data over all 2D dissociation constants
% (i.e. different ligands) and ligand copy-numbers (i.e. different doses)
y_min_Chmielewski = min(min(response_data_Chmielewski));
% Calculate the maximum response data over all 2D dissociation constants
% (i.e. different ligands) and ligand copy-numbers (i.e. different doses)
y_max_Chmielewski = max(max(response_data_Chmielewski));
% Calculate the extension to add to the y-axis
y_ext_Chmielewski = ax_scale*(y_max_Chmielewski-y_min_Chmielewski);
% Calculate the lower y-axis value
y_lower_Chmielewski = y_min_Chmielewski - y_ext_Chmielewski;
% Calculate the upper y-axis value
y_upper_Chmielewski = y_max_Chmielewski + y_ext_Chmielewski;

% Set a storage matrix for analysing the response data
response_Chmielewski = zeros(n_dose_Chmielewski, n_ligands_Chmielewski);
% Get response data from first ligand
response_Chmielewski(:,1) = data_Chmielewski.(char(data_Chmielewski.Properties.VariableNames(2)));
% Calculate a vector of NaN values where 0 represents NaN and 1 otherwise
index_nan_Chmielewski = ~isnan(response_Chmielewski(:,1));
% Set the second subplot
subplot(2,2,2)
% Plot dose against response for first ligand (for all non NaN values)
plot(dose_Chmielewski(index_nan_Chmielewski), response_Chmielewski(index_nan_Chmielewski,1),...
     'Marker', '*', 'Color', plot_col(1,:));
 % Loop through remaining ligands and plot dose response as above
for i=1:(n_ligands_Chmielewski-1)
    response_Chmielewski(:,i+1) = data_Chmielewski.(char(data_Chmielewski.Properties.VariableNames(i+2)));
    index_nan_Chmielewski = ~isnan(response_Chmielewski(:,i+1));
    line(dose_Chmielewski(index_nan_Chmielewski), response_Chmielewski(index_nan_Chmielewski,i+1),...
        'Marker', '*', 'Color', plot_col(i+1,:))
end
% Set the x-axis limits
xlim([x_lower_Chmielewski, x_upper_Chmielewski])
% Set the y-axis limits
ylim([y_lower_Chmielewski,y_upper_Chmielewski])
% Add the panel letter
text(x_lower_Chmielewski + pan_let_x_scale*(x_upper_Chmielewski-x_lower_Chmielewski),...
    y_lower_Chmielewski + pan_let_y_scale*(y_upper_Chmielewski-y_lower_Chmielewski),...
    '\fontsize{12} \bf b')
% Get the current axis handle
ax = gca;
% Set the position of the axes
ax.Position = [plot_sub_x_right, plot_sub_y_top, plot_sub_width, plot_sub_height];
% Add the legend text
l = legend(leg_txt_Chmielewski, 'Interpreter', 'latex', 'Location', 'best');
% Add the legend title
title(l, leg_tit, 'Interpreter', 'latex')
% Remove the box surrounding the legend
legend('boxoff')
% Add the plot title
title('Chmielewski $\emph{et al.}$ 2004', 'Interpreter', 'latex')
% Add the x-axis label
xlabel(x_lab, 'Interpreter', 'latex')
% Add the y-axis label
ylabel('IFN$\gamma$ (ng/ml)', 'Interpreter', 'latex')

%%% Plot mcmahan2006 data

% Set the 3D Kd values taken from mcmahan2006
K_d_3D_McMahan = [1.2, 1.9, 2.7, 3.1, 3.8, 4.1, 6.6];
% Set the legend text
leg_txt_McMahan = {num2str(round(log10(K_d_3D_McMahan),1)')};

% Set the name of the file containing the data
file_name_McMahan = 'data_mcmahan_2006.csv';
% Read the data
data_McMahan = readtable(file_name_McMahan, 'ReadVariableNames', 1);

% Calculate the total number of different ligands
n_ligands_McMahan = length(data_McMahan.Properties.VariableNames) - 1;
% Calculate dose values of a log10 scale
dose_McMahan = log10(data_McMahan.(char(data_McMahan.Properties.VariableNames(1))));
% Calculate the total number of dose values
n_dose_McMahan = length(dose_McMahan);

% Calculate the maximum dose
x_max_McMahan = max(dose_McMahan);
% Calculate the minimum dose
x_min_McMahan = min(dose_McMahan);
% Calculate the extension to add to the x-axis
x_ext_McMahan = ax_scale*(x_max_McMahan-x_min_McMahan);
% Calculate the lower x-axis value
x_lower_McMahan = -10; % Hard-coded to fit legend in
% Calculate the upper x-axis value
x_upper_McMahan = x_max_McMahan + x_ext_McMahan;

% Get all of the response data
% NOTE THE CONVERSION FROM per min to per sec
response_data_McMahan = permin2persec*table2array(data_McMahan(:,2:end));
% Calculate the minimum response data over all 2D dissociation constants
% (i.e. different ligands) and ligand copy-numbers (i.e. different doses)
y_min_McMahan = min(min(response_data_McMahan));
% Calculate the maximum response data over all 2D dissociation constants
% (i.e. different ligands) and ligand copy-numbers (i.e. different doses)
y_max_McMahan = max(max(response_data_McMahan));
% Calculate the extension to add to the y-axis
y_ext_McMahan = ax_scale*(y_max_McMahan-y_min_McMahan);
% Calculate the lower y-axis value
y_lower_McMahan = y_min_McMahan - y_ext_McMahan;
% Calculate the upper y-axis value
y_upper_McMahan = y_max_McMahan + y_ext_McMahan;

% Set a storage matrix for analysing the response data
response_McMahan = zeros(n_dose_McMahan, n_ligands_McMahan);
% Get response data from first ligand
% NOTE THE CONVERSION FROM per min to per sec
response_McMahan(:,1) = permin2persec*data_McMahan.(char(data_McMahan.Properties.VariableNames(2)));
% Calculate a vector of NaN values where 0 represents NaN and 1 otherwise
index_nan_McMahan = ~isnan(response_McMahan(:,1));
% Set the third subplot
subplot(2,2,3)
% Plot dose against response for first ligand (for all non NaN values)
plot(dose_McMahan(index_nan_McMahan), response_McMahan(index_nan_McMahan,1),...
     'Marker', '*', 'Color', plot_col(1,:));
% Loop through remaining ligands and plot dose response as above
for i=1:(n_ligands_McMahan-1)
    % NOTE THE CONVERSION FROM per min to per sec
    response_McMahan(:,i+1) = permin2persec*data_McMahan.(char(data_McMahan.Properties.VariableNames(i+2)));
    index_nan_McMahan = ~isnan(response_McMahan(:,i+1));
    line(dose_McMahan(index_nan_McMahan), response_McMahan(index_nan_McMahan,i+1),...
       'Marker', '*', 'Color', plot_col(i+1,:))
end
% Set the x-axis limits
xlim([x_lower_McMahan, x_upper_McMahan])
% Set the y-axis limits
ylim([y_lower_McMahan,y_upper_McMahan])
% Add the panel letter
text(x_lower_McMahan + pan_let_x_scale*(x_upper_McMahan-x_lower_McMahan),...
    y_lower_McMahan + pan_let_y_scale*(y_upper_McMahan-y_lower_McMahan),...
    '\fontsize{12} \bf c')
% Get the current axis handle
ax = gca;
% Set the position of the axes
ax.Position = [plot_sub_x_left, plot_sub_y_bottom, plot_sub_width, plot_sub_height];
% Add the legend text
l = legend(leg_txt_McMahan, 'Interpreter', 'latex', 'Location', 'best');
% Add the legend title
title(l, leg_tit, 'Interpreter', 'latex')
% Remove the box surrounding the legend
legend('boxoff')
% Add the plot title
title('McMahan $\emph{et al.}$ 2006', 'Interpreter', 'latex')
% Add the x-axis label
xlabel(x_lab, 'Interpreter', 'latex')
% Add the y-axis label
ylabel('T cell proliferation (cycles/s)', 'Interpreter', 'latex')

%%% Plot mean of previously unpublished Abu-Shah data

% Set the 3D Kd values taken from Pettmann2021
K_d_3D = [4.93, 22.8, 47.3, 162, 140, 299];
% Set the legend text
leg_txt_AbuShah = {num2str(round(log10(K_d_3D),1)')};

% Set the name of the file containing the data
file_name_AbuShah = 'data_abushah_2021.xlsx';
% Set the number of repeated experiments
n_reps_AbuShah = 3;
% Set the names of the sheets in the file
sheet1_AbuShah = 'rep1';
sheet2_AbuShah = 'rep2';
sheet3_AbuShah = 'rep3';
% Read the data from the different sheets
data1_AbuShah = readtable(file_name_AbuShah, 'Sheet', sheet1_AbuShah,...
    'ReadVariableNames', 1);
data2_AbuShah = readtable(file_name_AbuShah, 'Sheet', sheet2_AbuShah,...
    'ReadVariableNames', 1);
data3_AbuShah = readtable(file_name_AbuShah, 'Sheet', sheet3_AbuShah,...
    'ReadVariableNames', 1);
% Calculate a vector that contains the number of rows and columns of
% the dose response data
size_data1_AbuShah = size(data1_AbuShah);
% Set a storage array for holding all of the dose-response data across the
% experimetal repeats
data_AbuShah_tmp = zeros(size_data1_AbuShah(1),size_data1_AbuShah(2),n_reps_AbuShah);
% Allocate each dose-response data-set to the storage array
data_AbuShah_tmp(:,:,1) = table2array(data1_AbuShah);
data_AbuShah_tmp(:,:,2) = table2array(data2_AbuShah);
data_AbuShah_tmp(:,:,3) = table2array(data3_AbuShah);
% Calculate the mean response across the experimental repeats
data_AbuShah = array2table(mean(data_AbuShah_tmp,n_reps_AbuShah,'omitnan'),...
    'VariableNames', data1_AbuShah.Properties.VariableNames);

% Calculate the number of different ligands
n_ligands_AbuShah = length(data_AbuShah.Properties.VariableNames) - 1;
% Calculate the dose values on a log10 scale
dose_AbuShah = log10(data_AbuShah.(char(data_AbuShah.Properties.VariableNames(1))));
% Calculate the total number of dose values
n_dose_AbuShah = length(dose_AbuShah);

% Calculate the maximum dose
x_max_AbuShah = max(dose_AbuShah);
% Calculate the minimum dose
x_min_AbuShah = min(dose_AbuShah);
% Calculate the extension to add to the x-axis
x_ext_AbuShah = ax_scale*(x_max_AbuShah-x_min_AbuShah);
% Calculate the lower x-axis value
x_lower_AbuShah = -5.25; % Hard-coded to fit in the legend
% Calculate the upper x-axis value
x_upper_AbuShah = x_max_AbuShah + x_ext_AbuShah;

% Get all of the response data 
% NOTE THE CONVERSION FROM pg/ml to ng/ml
response_data_AbuShah = pg2ng*table2array(data_AbuShah(:,2:end));
% Calculate the minimum response data over all 2D dissociation constants
% (i.e. different ligands) and ligand copy-numbers (i.e. different doses)
y_min_AbuShah = min(min(response_data_AbuShah));
% Calculate the maximum response data over all 2D dissociation constants
% (i.e. different ligands) and ligand copy-numbers (i.e. different doses)
y_max_AbuShah = max(max(response_data_AbuShah));
% Calculate the extension to add to the y-axis
y_ext_AbuShah = ax_scale*(y_max_AbuShah-y_min_AbuShah);
% Calculate the lower y-axis value
y_lower_AbuShah = y_min_AbuShah - y_ext_AbuShah;
% Calculate the upper y-axis value
y_upper_AbuShah = y_max_AbuShah + y_ext_AbuShah;

% Set a storage matrix for analysing the response data
response_AbuShah = zeros(n_dose_AbuShah, n_ligands_AbuShah);
% Get response data from first ligand
% NOTE THE CONVERSION FROM pg/ml to ng/ml
response_AbuShah(:,1) = pg2ng*data_AbuShah.(char(data_AbuShah.Properties.VariableNames(2)));
% Calculate a vector of NaN values where 0 represents NaN and 1 otherwise
index_nan_AbuShah = ~isnan(response_AbuShah(:,1));
% Set the fourth subplot
subplot(2,2,4)
% Plot dose against response for first ligand (for all non NaN values)
plot(dose_AbuShah(index_nan_AbuShah), response_AbuShah(index_nan_AbuShah,1),...
     'Marker', '*', 'Color', plot_col(1,:));
% Loop through remaining ligands and plot dose response as above
for i=1:(n_ligands_AbuShah-1)
    % NOTE THE CONVERSION FROM pg/ml to ng/ml
    response_AbuShah(:,i+1) = pg2ng*data_AbuShah.(char(data_AbuShah.Properties.VariableNames(i+2)));
    index_nan_AbuShah = ~isnan(response_AbuShah(:,i+1));
    line(dose_AbuShah(index_nan_AbuShah), response_AbuShah(index_nan_AbuShah,i+1),...
       'Marker', '*', 'Color', plot_col(i+1,:))
end
% Set the x-axis limits 
xlim([x_lower_AbuShah, x_upper_AbuShah])
% Set the y-axis limits
ylim([y_lower_AbuShah,y_upper_AbuShah])
% Add the panel letter
text(x_lower_AbuShah + pan_let_x_scale*(x_upper_AbuShah-x_lower_AbuShah),...
    y_lower_AbuShah + pan_let_y_scale*(y_upper_AbuShah-y_lower_AbuShah),...
    '\fontsize{12} \bf d')
% Get the current axis handle
ax = gca;
% Set the position of the axes
ax.Position = [plot_sub_x_right, plot_sub_y_bottom, plot_sub_width, plot_sub_height];
% Set the position of the axes
l = legend(leg_txt_AbuShah, 'Interpreter', 'latex', 'Location', 'best');
% Add the plot title
title(l, leg_tit, 'Interpreter', 'latex')
% Remove the box surrounding the legend
legend('boxoff')
% Add the plot title
title('Previously unpublished data', 'Interpreter', 'latex')
% Add the x-axis label
xlabel(x_lab, 'Interpreter', 'latex')
% Add the y-axis label
ylabel('IL-2 (ng/ml)', 'Interpreter', 'latex')

%%% Write figure to file

% Save as a .eps file
print('fig_dose_resp_expr', '-depsc')
% Save as a .emf file
print('fig_dose_resp_expr', '-dmeta')
