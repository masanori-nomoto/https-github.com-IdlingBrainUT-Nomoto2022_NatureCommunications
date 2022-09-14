%% figure line plot
% [ttest_results] = fxn_line_plot(data1_all, data2_all, xtickname1, xtickname2, ylabel_title_input) 
clc; clear; close all

load('NatComm_data_Fig6ij_mahal_ca1_ca3_cs_iti_us_iti_all.mat')
%% data input

% ### for cs
data1_all = mahal_ca3_cs_iti_all_wt;
data2_all = mahal_ca3_cs_iti_all_ko;

xtickname1 = 'CS'
xtickname2 = 'ITI'
ylabel_title_input =  {'Relative change in';'PVD to ITI (%)'};
%%
[ttest_results] = fxn_line_plot_for_imported_data(data1_all, data2_all, xtickname1, xtickname2, ylabel_title_input) ;
%%

% ### for us
data1_all = mahal_ca3_us_iti_all_wt;
data2_all = mahal_ca3_us_iti_all_ko;

xtickname1 = 'US'
xtickname2 = 'ITI'
ylabel_title_input =  {'Relative change in';'PVD to ITI (%)'};

close all
%%
[ttest_results] = fxn_line_plot_for_imported_data(data1_all, data2_all, xtickname1, xtickname2, ylabel_title_input) ;
%%