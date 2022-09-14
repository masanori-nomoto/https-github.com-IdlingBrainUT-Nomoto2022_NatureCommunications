%% CA3 triple connectivity summarized code 211021version
% for tiple
%% load mat-file
clc; clear; close all;
load('NatComm_data_CA3_calcium_data_raw_ctrl11_ko9.mat');
tic;
%% #### select this session or below #### 
% data loading c
input_data = {};
input_data_num = 11;

input_data{1,1} = ca3_ctrl_m01;
input_data{2,1} = ca3_ctrl_m02;
input_data{3,1} = ca3_ctrl_m03;
input_data{4,1} = ca3_ctrl_m04;
input_data{5,1} = ca3_ctrl_m05;
input_data{6,1} = ca3_ctrl_m06;
input_data{7,1} = ca3_ctrl_m07;
input_data{8,1} = ca3_ctrl_m08;
input_data{9,1} = ca3_ctrl_m09;
input_data{10,1} = ca3_ctrl_m10;
input_data{11,1} = ca3_ctrl_m11;

disp('finish data loading')
%% #### select this session or below #### 
input_data = {};
input_data_num = 9;

input_data{1,1} = ca3_ko_m01;
input_data{2,1} = ca3_ko_m02;
input_data{3,1} = ca3_ko_m03;
input_data{4,1} = ca3_ko_m04;
input_data{5,1} = ca3_ko_m05;
input_data{6,1} = ca3_ko_m06;
input_data{7,1} = ca3_ko_m07;
input_data{8,1} = ca3_ko_m08;
input_data{9,1} = ca3_ko_m09;

disp('finish data loading')
%% calculate results.out
% 10: shared, 11: non-shared
% sigma_cut_if_ON = 3; % dafault, Z score cutoff 
% [result]  = fxn_CA3_connectivity_for_summary(input_data)
out_table = {};
for i = 1:input_data_num
data_temp = cell2mat(input_data(i,1));

% ### for all session
out_table{i,1}  = fxn_CA3_triple_connectivity_cal_summary(data_temp); % for triple, select here or below by comment out

end
toc;
%%
 bin_summary = {};  bin_ITI1 = {}; bin_ITI2 = {};
 bin_500ms_summary  = {}; bin_500ms_ITI1  = {}; bin_500ms_ITI2 = {};
 bin_500ms_summary{10,2} = 'shared'; bin_500ms_summary{11,2} = 'non-shared';
for i_big = 1:11 % for triple, 10: shared, 11: non-shared in triple
    for i = 1:input_data_num
%    bin_summary{i_big,1}(i,:)        = out_table{i,1}.summary(i_big,:);
%    bin_ITI1{i_big,1}(i,:)            = out_table{i,1}.ITI1(i_big,:);
%    bin_ITI2{i_big,1}(i,:)             = out_table{i,1}.ITI2(i_big,:);
   bin_500ms_summary{i_big,1}(i,:)    = out_table{i,1}.summary_500ms(i_big,:);      %  for triple-wise synchrony
   bin_500ms_sync{i_big,1}(i,:)          = out_table{i, 1}.total_sync(i_big,:);     %  for event rate (hz)
   bin_500ms_cell{i_big,1}(i,:)      = out_table{i, 1}.total_cell_500ms(i_big,:);   %  for maximal combination of cells
   bin_500ms_time{i_big,1}(i,:)      = out_table{i, 1}.time_500ms(i_big,:);         %  for divided frames in each session
   bin_500ms_triple_sync_raw_500ms{i_big,1}(i,:)      = out_table{i, 1}.triple_sync_raw_500ms(i_big,:);  % for actual synchrony count
   bin_500ms_triple_sync_standard_500ms{i_big,1}(i,:)      = out_table{i, 1}.triple_sync_standard_500ms(i_big,:); % for normalized synchrony count by max combination of cells
%    bin_500ms_ITI1{i_big,1}(i,:)      = out_table{i,1}.ITI1_500ms(i_big,:);
%    bin_500ms_ITI2{i_big,1}(i,:)     = out_table{i,1}.ITI2_500ms(i_big,:);
%    bin_500ms_summary_shuffled{i_big,1}(i,:) = out_table{i,1}.summary_500ms_shuffled(i_big,:);
    end
end
toc;
%% Results for figure 6e, 6f, and 6g

% bin_500ms_summary is for synchrony index: shared (10th raw) and non-shared (11th raw) ;
% bin_500ms_sync    is for synchrony rate : non-shared (11th raw);

% column; 1st: baseline, 2nd: CS, 3rd: US, 
% 4th: ITI-1, 5th: ITI-2, 6th: rest.

%%