%% Boxplot code for import data
%% for final 500ms data
clc; clear; close all;

% ### for iti1-15s data
load('NatComm_data_Fig3ef_500ms_iti1_15s_shared_non_shared_connectivity.mat')
input_cell = ca1_revision_500ms_iti1_15s_shared_non_shared_connectivity;

ctrl_range = [1:6]; ko_range = [10:20]; % CA3 data, shared, and non-shared
%% figure boxplot for shared subpopulation
% fxn_boxplot_double(input_cell, ctrl_range, ko_range, yaxis_max);
% for shared
% yaxis_max = 0.01; % for log scale
% yaxis_max = 0.1; % for log version
yaxis_max = 0.01; % 
comparison_section = 1; disp('Connectivity between shared analyzed')% CS vs LTM-CS

fxn_boxplot_double_for_ca3_revised_importdata(input_cell, yaxis_max, comparison_section, ctrl_range, ko_range); box off
% set(gca,'YScale','log'); % for ca3 triple shared, log
ylim([0 0.002]);  %yticks([0 0.001 0.002])
% ax = gca;
% ax.YAxis.Exponent = 5;
%% figure boxplot for non-shared subpopulation
% fxn_boxplot_double(input_cell, ctrl_range, ko_range, yaxis_max);
% for non-shared
% yaxis_max = 0.01; % for log scale
yaxis_max = 0.0001; % for log version
comparison_section = 2; disp('Connectivity between non-shared analyzed')% CS vs LTM-CS

fxn_boxplot_double_for_ca3_revised_importdata(input_cell, yaxis_max, comparison_section, ctrl_range, ko_range); box off
ylim([0 0.0002])
 %%