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
out_table{i,1}  = fxn_CA1_connectivity_cal_summary(data_temp);
end
toc;
%%
 bin_summary = {};  bin_ITI1 = {}; bin_ITI2 = {};
 bin_500ms_summary  = {}; bin_500ms_ITI1  = {}; bin_500ms_ITI2 = {};
 bin_500ms_summary{10,2} = ('Supplementary Fig4c');
for i_big = 1:10 % 10th row is AB vs. pure-C data!
    for i = 1:input_data_num
   bin_summary{i_big,1}(i,:)        = out_table{i,1}.summary(i_big,:);
   bin_ITI1{i_big,1}(i,:)            = out_table{i,1}.ITI1(i_big,:);
   bin_ITI2{i_big,1}(i,:)             = out_table{i,1}.ITI2(i_big,:);
   bin_500ms_summary{i_big,1}(i,:)    = out_table{i,1}.summary_500ms(i_big,:); % 10th row is AB vs. pure-C data!
%    bin_500ms_ITI1{i_big,1}(i,:)      = out_table{i,1}.ITI1_500ms(i_big,:);
%    bin_500ms_ITI2{i_big,1}(i,:)     = out_table{i,1}.ITI2_500ms(i_big,:);

%    bin_500ms_sync{i_big,1}(i,:)      = out_table{i, 1}.total_sync(i_big,:);
%    bin_500ms_cell{i_big,1}(i,:)      = out_table{i, 1}.total_cell_500ms(i_big,:);
%    bin_500ms_time{i_big,1}(i,:)      = out_table{i, 1}.time_500ms(i_big,:);   
    end
end
toc;
%% #### for Supplementary Fig4c ####

% bin_500ms_summary: % 10th row is ITI-1 synchrony between CS-US vs. pure Test-CS cells data!

%%