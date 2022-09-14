%% CA1 connectivity summarized code 
%% load mat-file
clc; clear; close all;
load('NatComm_data_CA1_calcium_data_raw_ctrl6_ko11.mat');
tic;
%% #### select this session or below #### 
% data loading c
input_data = {};
input_data_num = 6;

input_data{1,1} = ca1_ctrl_m01;
input_data{2,1} = ca1_ctrl_m02;
input_data{3,1} = ca1_ctrl_m03;
input_data{4,1} = ca1_ctrl_m04;
input_data{5,1} = ca1_ctrl_m05;
input_data{6,1} = ca1_ctrl_m06;

disp('finish data loading')
%% #### select this session or below #### 
input_data = {};
input_data_num = 11;

input_data{1,1} = ca1_ko_m01;
input_data{2,1} = ca1_ko_m02;
input_data{3,1} = ca1_ko_m03;
input_data{4,1} = ca1_ko_m04;
input_data{5,1} = ca1_ko_m05;
input_data{6,1} = ca1_ko_m06;
input_data{7,1} = ca1_ko_m07;
input_data{8,1} = ca1_ko_m08;
input_data{9,1} = ca1_ko_m09;
input_data{10,1} = ca1_ko_m10;
input_data{11,1} = ca1_ko_m11;

disp('finish data loading')

%% calculate results.out
% [result]  = fxn_CA1_connectivity_for_summary(input_data)
out_table = {};
for i = 1:input_data_num
data_temp = cell2mat(input_data(i,1));
out_table{i,1}  = fxn_CA1_triple_connectivity_cal_summary(data_temp); % normal
% out_table{i,1}  = fxn_CA1_triple_connectivity_cal_summary_shuffle(data_temp); % for shuffling
end
toc;
%%
% for triple, 10: shared, 11: non-shared in triple
 bin_summary = {};  bin_ITI1 = {}; bin_ITI2 = {};
 bin_500ms_summary  = {}; bin_500ms_ITI1  = {}; bin_500ms_ITI2 = {};
 bin_500ms_summary{10,2} = 'shared'; bin_500ms_summary{11,2} = 'non-shared';

for i_big = 1:11 % for triple, 10: shared, 11: not-shared in triple
    for i = 1:input_data_num
%    bin_summary{i_big,1}(i,:)                = out_table{i,1}.summary(i_big,:);
%    bin_ITI1{i_big,1}(i,:)                 = out_table{i,1}.ITI1(i_big,:);
%    bin_ITI2{i_big,1}(i,:)                 = out_table{i,1}.ITI2(i_big,:);
   bin_500ms_summary{i_big,1}(i,:)          = out_table{i,1}.summary_500ms(i_big,:);
% 220818 revision
   bin_500ms_sync{i_big,1}(i,:)          = out_table{i, 1}.total_sync(i_big,:);
   bin_500ms_cell{i_big,1}(i,:)      = out_table{i, 1}.total_cell_500ms(i_big,:);
   bin_500ms_time{i_big,1}(i,:)      = out_table{i, 1}.time_500ms(i_big,:);
   bin_500ms_triple_sync_raw_500ms{i_big,1}(i,:)      = out_table{i, 1}.triple_sync_raw_500ms(i_big,:);
   bin_500ms_triple_sync_standard_500ms{i_big,1}(i,:)      = out_table{i, 1}.triple_sync_standard_500ms(i_big,:);
%    bin_500ms_ITI1{i_big,1}(i,:)           = out_table{i,1}.ITI1_500ms(i_big,:);
%    bin_500ms_ITI2{i_big,1}(i,:)           = out_table{i,1}.ITI2_500ms(i_big,:);
%    bin_500ms_summary_shuffled{i_big,1}(i,:) = out_table{i,1}.summary_500ms_shuffled(i_big,:);
    end
end
toc;
%% ### Results for figure 3e, 3f, and 6g (for CA1 ITI-1 synchrony rate raw count) ###

% bin_500ms_summary is for synchrony index: shared (10th raw) and non-shared (11th raw) ;
% bin_500ms_sync    is for synchrony rate : non-shared (11th raw);

% column; 1st: baseline, 2nd: CS, 3rd: US, 
% 4th: ITI-1, 5th: ITI-2, 6th: rest.

% bin_500ms_summary_supple4c is for supplementary fig 4c.

%%



