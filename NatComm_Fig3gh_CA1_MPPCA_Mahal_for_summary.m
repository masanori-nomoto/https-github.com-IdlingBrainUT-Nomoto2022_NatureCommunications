%% CA1 MPPCA-Mahal summarized code 
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

%% CA1 MP-PCA Mahalanobis distance analysis (210609 revised)
%  [res_MPPCA_data_cell, res_Angle_data_cell, res_Mahald_data_cell]  = ...
%     fxn_CA1_MP_PCA_Mahalanobis_summary(input_data, thrcov_PC_percent, reference_session_num, num_top_PCs)
thrcov_PC_percent = 70 ; % 
analyzed_session_num = 5 ; % 2:base, 3:CS, 4:US, 5:ITI1, 6:ITI2, 7:LTM_base, 8:LTM-CS, 9:LTM-ITI
num_top_PCs = 3 ; % 3 off

out_table = {};
for i = 1:input_data_num
data_temp = cell2mat(input_data(i,1));
[res_MPPCA_data_cell{i,1}, res_Angle_data_cell{i,1}, res_Mahald_data_cell{i,1}, res_CS_cells_data_cell{i,1}] = ...
    fxn_CA1_MP_PCA_Mahalanobis_summary(data_temp, thrcov_PC_percent, analyzed_session_num, num_top_PCs);
end
% out_table{i,1}  = 
toc;
%%
%  bin_summary = {};  bin_ITI1 = {}; bin_ITI2 = {};
%  bin_200ms_summary  = {}; bin_200ms_ITI1  = {}; bin_200ms_ITI2 = {};
Mahald_mean_AB  = {}; Mahald_mean_AC  = {}; Mahald_mean_BC  = {}; Mahald_mean_ABC = {};
Mahald_mean_A = {}; Mahald_mean_B = {}; Mahald_mean_C = {};

for i = 1:input_data_num
Mahald_mean{i,1} = res_Mahald_data_cell{i, 1}(11:end,2:end);
end

for i = 1:input_data_num
% Mahald_mean_AB{i,1}  = cell2mat(Mahald_mean{i,1}(1,:));
% Mahald_mean_AC{i,1}  = cell2mat(Mahald_mean{i,1}(2,:));
% Mahald_mean_BC{i,1}  = cell2mat(Mahald_mean{i,1}(3,:));
% Mahald_mean_ABC{i,1} = cell2mat(Mahald_mean{i,1}(4,:));
Mahald_mean_A{i,1}  = cell2mat(Mahald_mean{i,1}(5,:));
Mahald_mean_B{i,1}  = cell2mat(Mahald_mean{i,1}(6,:));
% Mahald_mean_C{i,1} = cell2mat(Mahald_mean{i,1}(7,:));
end

% Mahald_mean_AB = cell2mat(Mahald_mean_AB);
% Mahald_mean_AC = cell2mat(Mahald_mean_AC);
% Mahald_mean_BC = cell2mat(Mahald_mean_BC);
% Mahald_mean_ABC = cell2mat(Mahald_mean_ABC);
Mahald_mean_A = cell2mat(Mahald_mean_A);
Mahald_mean_B = cell2mat(Mahald_mean_B);
% Mahald_mean_C = cell2mat(Mahald_mean_C);
%%
%  bin_summary = {};  bin_ITI1 = {}; bin_ITI2 = {};
%  bin_200ms_summary  = {}; bin_200ms_ITI1  = {}; bin_200ms_ITI2 = {};
Angle_mean_AB  = {};
Angle_mean_AC  = {};
Angle_mean_BC  = {};
Angle_mean_ABC = {};
Angle_mean_A = {}; Angle_mean_B = {}; Angle_mean_C = {};

for i = 1:input_data_num
Angle_mean{i,1} = res_Angle_data_cell{i, 1}(11:end,2:end);
end

for i = 1:input_data_num
% Angle_mean_AB{i,1}  = cell2mat(Angle_mean{i,1}(1,:));
% Angle_mean_AC{i,1}  = cell2mat(Angle_mean{i,1}(2,:));
% Angle_mean_BC{i,1}  = cell2mat(Angle_mean{i,1}(3,:));
% Angle_mean_ABC{i,1} = cell2mat(Angle_mean{i,1}(4,:));
Angle_mean_A{i,1}  = cell2mat(Angle_mean{i,1}(5,:));
Angle_mean_B{i,1}  = cell2mat(Angle_mean{i,1}(6,:));
% Angle_mean_C{i,1} = cell2mat(Angle_mean{i,1}(7,:));
end

% Angle_mean_AB = cell2mat(Angle_mean_AB);
% Angle_mean_AC = cell2mat(Angle_mean_AC);
% Angle_mean_BC = cell2mat(Angle_mean_BC);
% Angle_mean_ABC = cell2mat(Angle_mean_ABC);
Angle_mean_A = cell2mat(Angle_mean_A);
Angle_mean_B = cell2mat(Angle_mean_B);
% Angle_mean_C = cell2mat(Angle_mean_C);
%%
for i = 1:input_data_num
CS_cells_data_mean(i,:) = cell2mat(res_CS_cells_data_cell{i, 1}(1,1));
end
%%
%% ### for figure, Relative PVD cal between CS or US and ITI-1 ###
CS_PVD_temp           = Mahald_mean_A(:,[4,2]);
CS_PVD_temp           = CS_PVD_temp ./ CS_PVD_temp(:,1);
results_PVD_CS_ITI_percent = (CS_PVD_temp-1)*100; % for fig. 3g

US_PVD_temp           = Mahald_mean_B(:,[4,3]);
US_PVD_temp           = US_PVD_temp ./ US_PVD_temp(:,1);
results_PVD_US_ITI_percent = (US_PVD_temp-1)*100; % for fig. 3h

%% ### for figure, Relative rotation cal between CS or US and ITI-1 ###
results_PVD_CS_ITI_angle = Angle_mean_A(:,[4,2]); % for fig3g CS angle
results_PVD_US_ITI_angle = Angle_mean_B(:,[4,3]); % for fig3h US angle

%%
