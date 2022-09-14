%% figure line plot
% [ttest_results] = fxn_line_plot(data1_all, data2_all, xtickname1, xtickname2, ylabel_title_input) 
close all; clc; clear;

load('NatComm_data_Fig3ef_for_Pearson_correlation_fig_freeze_connectivity_index_in_ca1.mat');

%% input data for non-shared selection

% ### data input triple not-shared
data1_all = triple_notshared_ctrl;
data2_all = triple_notshared_ko;

%%
data1_all = [data1_all(:,1).*100, data1_all(:,2)]; % conrevisedCA3500ms1sconnectivityITI2nectivity % conversion
data2_all = [data2_all(:,1).*100, data2_all(:,2)]; % connectivity % conversion

xlabel_title_input =  {'Change in ITI connectivity';'to baseline (%)'};
ylabel_title_input =  {'Freezing (%)'};
% ylabel_title_input =  {'Freezing';'in CS (%)'};
%%
fig_xy_dim = [140 125]; % input [x y] dim in pixel
fxn_show_scatter_mu_error_linearreg(fig_xy_dim, data1_all, data2_all, xlabel_title_input, ylabel_title_input) ;
%%

%% input data for shared selection

% ### data input CS-US
data1_all = triple_shared_ctrl;
data2_all = triple_shared_ko;
%%
data1_all = [data1_all(:,1).*100, data1_all(:,2)]; % conrevisedCA3500ms1sconnectivityITI2nectivity % conversion
data2_all = [data2_all(:,1).*100, data2_all(:,2)]; % connectivity % conversion

xlabel_title_input =  {'Change in ITI connectivity';'to baseline (%)'};
ylabel_title_input =  {'Freezing (%)'};
% ylabel_title_input =  {'Freezing';'in CS (%)'};
%%
fig_xy_dim = [140 125]; % input [x y] dim in pixel
fxn_show_scatter_mu_error_linearreg(fig_xy_dim, data1_all, data2_all, xlabel_title_input, ylabel_title_input) ;
%%