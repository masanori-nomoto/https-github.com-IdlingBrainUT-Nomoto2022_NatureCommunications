%% CA1 MP-PCA Mahalanobis distance analysis (210609 revised)
function [res_MPPCA_data_cell, res_Angle_data_cell, res_Mahald_data_cell, res_CS_cells_data_cell]  = ...
    fxn_CA1_MP_PCA_Mahalanobis_summary(input_data, thrcov_PC_percent, reference_session_num, num_top_PCs)
% clc; clear; close all; disp('Previous data are cleared!'); 
%% for debug,
% thrcov_PC_percent = 50;
% reference_session_num = 5; % 5: ITI1, positive results!
%%
disp('MP-PCA-Mahalanobis cal starts!');
%% determine connectivity calculating rate
% upsample_rate = 1; % for 1s binning data
% upsample_rate = 4; % for 200ms binning data, 4 x 50ms = 200ms.

% for 1s binning data
ca_binning_num = 20; % 20Hz to 1Hz data
timestamp_upsample_rate = 1; % upsaple with factor 1 (no upsample)

% for 200ms data
% ca_binning_num = 4; % 20Hz to 1Hz data
% timestamp_upsample_rate = 4;

%% filtering ans sigma-cut off
% ca_filt_data = fxn_filt_sigma_ca(ca_raw_data, sigma_if_z_ON);

sigma_cut_if_ON = 3; % Z score cutoff 
G1_concatenated_temp  = fxn_filt_sigma_cutoff(input_data, sigma_cut_if_ON);
%% Temporal binning
% [Ca_data_binned,  Ca_data_binned_mean] = funcHF_temporal_binning(Ca_data, binning_num)
binning_num = 20; % 20Hz to 1Hz data
[G1_Ca_binned,  G1_Ca_binned_mean] = fxn_temporal_binning(G1_concatenated_temp, binning_num);

reshape_range = (16:600);
% reshape_range = (151:600);
fxn_1s_bin_matrix(G1_Ca_binned, G1_Ca_binned, G1_Ca_binned_mean, G1_Ca_binned, reshape_range)  
%% input time stamp information
% upsample_rate = 1; % for 1s binning data
% upsample_rate = 4; % for 200ms binning data, 4 x 50ms = 200ms.

time_stamp_info = fxn_ca62400_12sessions(timestamp_upsample_rate); % load session time information
total_session_num  = 12;
responsiveness_threshold = 2;

G1_results = time_stamp_info;
%% Calculate mean in each session
for i = 1:total_session_num
G1_results{i,2} = G1_Ca_binned(time_stamp_info{i,1},:);
G1_results{i,3} = mean(G1_results{i,2},1);
G1_results{i,6} = G1_results{i,2} > 0;
end
%% Calculate responsiveness and cell-ID extraction
% find(X>0 & X<10,3)
for i = 2:total_session_num
G1_results{i,4} = (G1_results{i,3})./G1_results{1,3};
G1_results{i,5} = find(G1_results{i,4} >= responsiveness_threshold & isfinite(G1_results{i,4}));
G1_results{i,7} = find(isfinite(G1_results{i,4}));

    % previous code
    % for i = 2:total_session_num
    % G1_results{i,4} = (G1_results{i,3})./G1_results{1,3};
    % G1_results{i,5} = find(G1_results{i,4} >= responsiveness_threshold); % neuron ID
end
%% Synchronous activity 
% input parameter
connectivity_results_t = {};

connectivity_results_t{1,1} = ('CS vs LTM-CS'); 
connectivity_results_t{2,1} = ('US vs LTM-CS');
connectivity_results_t{3,1} = ('CS vs US');
connectivity_results_t{4,5} = ('ITI-1'); connectivity_results_t{4,6} = ('ITI-2');

for i_index_cal = 1:total_session_num;
data_input_temp = G1_results{i_index_cal,6};  % {2 is input_data_session (2=CS, 3=US, 4=ITI), fixed as 6}

Cell_type_A_ID  = G1_results{2,5};  % {2 is CS+ cell-IDs for refference, fixed as 5}
Cell_type_B_ID  = G1_results{3,5};  % {3 is US+ cell-IDs for refference, fixed as 5}
Cell_type_ITI1_ID  = G1_results{4,5};  % {3 is US+ cell-IDs for refference, fixed as 5}
Cell_type_ITI2_ID  = G1_results{5,5};  % {3 is US+ cell-IDs for refference, fixed as 5}
Cell_type_C_ID  = G1_results{11,5}; % {13 is LTM-CS+ cell-IDs for refference, fixed as 5} % homecage
% Cell_type_C_ID  = G1_results{13,5}; % {13 is LTM-CS+ cell-IDs for refference, fixed as 5}

    % 210524 added
    Cell_type_A_C_shared_ID   = intersect(Cell_type_A_ID,Cell_type_C_ID);
    Cell_type_A_C_not_shared_A_ID = setdiff(Cell_type_A_ID, Cell_type_A_C_shared_ID);
    Cell_type_A_C_not_shared_C_ID = setdiff(Cell_type_C_ID, Cell_type_A_C_shared_ID);

    Cell_type_A_B_shared_ID   = intersect(Cell_type_A_ID,Cell_type_B_ID);
    Cell_type_A_B_not_shared_A_ID = setdiff(Cell_type_A_ID, Cell_type_A_B_shared_ID);
    Cell_type_A_B_not_shared_B_ID = setdiff(Cell_type_B_ID, Cell_type_A_B_shared_ID);

    Cell_type_B_C_shared_ID   = intersect(Cell_type_B_ID,Cell_type_C_ID);
    Cell_type_B_C_not_shared_B_ID = setdiff(Cell_type_B_ID, Cell_type_B_C_shared_ID);
    Cell_type_B_C_not_shared_C_ID = setdiff(Cell_type_C_ID, Cell_type_B_C_shared_ID);
    
    Cell_type_A_B_C_shared_ID   = intersect(Cell_type_A_B_shared_ID,Cell_type_C_ID);

%% count number of shared-cell
Cell_merge_A_B = numel(intersect(Cell_type_A_ID,Cell_type_B_ID));
Cell_merge_A_C = numel(intersect(Cell_type_A_ID,Cell_type_C_ID));
Cell_merge_B_C = numel(intersect(Cell_type_B_ID,Cell_type_C_ID));
Cell_merbe_triple = numel(intersect(intersect(Cell_type_A_ID,Cell_type_B_ID),Cell_type_C_ID));

Cell_merge_ITI1_A = numel(intersect(Cell_type_A_ID,Cell_type_ITI1_ID));
Cell_merge_ITI1_C = numel(intersect(Cell_type_ITI1_ID,Cell_type_C_ID));
Cell_merbe_ITI_triple = numel(intersect(intersect(Cell_type_A_ID, Cell_type_ITI1_ID),Cell_type_C_ID));
%% Cell ID selection for MP-PCA cell resorting
MPPCA_data_cell = {};

MPPCA_data_cell{1,2} = 'Resort IDs'; MPPCA_data_cell{1,3} = 'Resort 1bin Ca'; 
MPPCA_data_cell{1,4} = 'MP-PCA'; MPPCA_data_cell{1,5} = 'thrcov_PCA';
MPPCA_data_cell{2,1} = 'AB'; MPPCA_data_cell{3,1} = 'AC'; MPPCA_data_cell{4,1} = 'BC'; MPPCA_data_cell{5,1} = 'ABC';
MPPCA_data_cell{6,1} = 'A'; MPPCA_data_cell{7,1} = 'B'; MPPCA_data_cell{8,1} = 'C';

MPPCA_data_cell{2,2} = unique([Cell_type_A_ID, Cell_type_B_ID]); % PCA_cellIDs_AB
MPPCA_data_cell{3,2} = unique([Cell_type_A_ID, Cell_type_C_ID]); % PCA_cellIDs_AC
MPPCA_data_cell{4,2} = unique([Cell_type_B_ID, Cell_type_C_ID]); % PCA_cellIDs_BC
MPPCA_data_cell{5,2} = unique([Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID]); % PCA_cellIDs_ABC
MPPCA_data_cell{6,2} = unique([Cell_type_A_ID]); % PCA_cellIDs_ABC
MPPCA_data_cell{7,2} = unique([Cell_type_B_ID]); % PCA_cellIDs_ABC
MPPCA_data_cell{8,2} = unique([Cell_type_C_ID]); % PCA_cellIDs_ABC
%% MP-PCA cal
%  [res_MPPCA, res_thrcov_PCA] = fxn_Marchenko_thrcover_PCA(ca_data, thrcov_PC_percent)
% thrcov_PC_percent = 90;
MPPCA_analysis_range = [1:600, 2761:3120];
% MPPCA_analysis_range = [1:600, 2401:2760, 2761:3120];

% for binary mode
% binary = find(G1_Ca_binned > 0);
% G1_Ca_binned(binary) = 1;

for i = 2:8
MPPCA_data_cell{i,3} = G1_Ca_binned(MPPCA_analysis_range,MPPCA_data_cell{i,2}) ; % PCA_cellIDs_AB
MPPCA_data_cell{i,3} = G1_Ca_binned(MPPCA_analysis_range,MPPCA_data_cell{i,2}) ; % PCA_cellIDs_AB
[MPPCA_data_cell{i,4}, MPPCA_data_cell{i,5}] = fxn_Marchenko_thrcover_PCA(MPPCA_data_cell{i,3}, thrcov_PC_percent); % G1_Ca_binned
end

%  [angle_raw, angle_mean] = fxn_nmt_cosine(reference_data, target_data)
%  [angle_raw, angle_mean] = fxn_nmt_cosine(reference_data, target_data)
%% Input cell info
Mahald_data_cell = {};

Mahald_data_cell{1,2} = 'Baseline'; Mahald_data_cell{1,3} = 'CS'; 
Mahald_data_cell{1,4} = 'US'; Mahald_data_cell{1,5} = 'ITI1'; Mahald_data_cell{1,6} = 'ITI2';
Mahald_data_cell{1,7} = 'LTM-Acc'; Mahald_data_cell{1,8} = 'LTM-CS'; Mahald_data_cell{1,9} = 'LTM-ITI';

Mahald_data_cell{2,1} = 'Time stamp';
Mahald_data_cell{2,2} = time_stamp_info{1,1}; Mahald_data_cell{2,3} = time_stamp_info{2,1};
Mahald_data_cell{2,4} = time_stamp_info{3,1}; Mahald_data_cell{2,5} = time_stamp_info{4,1}; 
Mahald_data_cell{2,6} = time_stamp_info{5,1};
Mahald_data_cell{2,7} = time_stamp_info{10,1}-2160;
Mahald_data_cell{2,8} = time_stamp_info{11,1}-2160;
Mahald_data_cell{2,9} = time_stamp_info{12,1}-2160;

Mahald_data_cell{3,1} = 'AB'; Mahald_data_cell{4,1} = 'AC'; Mahald_data_cell{5,1} = 'BC'; Mahald_data_cell{6,1} = 'ABC';
Mahald_data_cell{7,1} = 'A'; Mahald_data_cell{8,1} = 'B'; Mahald_data_cell{9,1} = 'C';
%% for cosine theta radian
Angle_data_cell = {};

Angle_data_cell{1,2} = 'Baseline'; Angle_data_cell{1,3} = 'CS'; 
Angle_data_cell{1,4} = 'US'; Angle_data_cell{1,5} = 'ITI1'; Angle_data_cell{1,6} = 'ITI2';
Angle_data_cell{1,7} = 'LTM-Acc'; Angle_data_cell{1,8} = 'LTM-CS'; Angle_data_cell{1,9} = 'LTM-ITI';

Angle_data_cell{2,1} = 'Time stamp';
Angle_data_cell{3,1} = 'AB'; Angle_data_cell{4,1} = 'AC'; Angle_data_cell{5,1} = 'BC'; Angle_data_cell{6,1} = 'ABC';
Angle_data_cell{7,1} = 'A'; Angle_data_cell{8,1} = 'B'; Angle_data_cell{9,1} = 'C';
%% Mahalanobis distance calculation
% [Mahalanobis_d_raw, Mahalanobis_mean] = fxn_nmt_Maharanobis(reference_data, target_data)
% MPPCA_data_cell{2, 5}.thrcov_PCA_score -> {5,5}

% reference_session = 2; % 2: baseline, 
% reference_session_num = 5; % 5: ITI1, positive results!
% reference_session = 8; % 8: LTM-CS, 

for ii = 2:8
for i = 2:9
[Mahald_data_cell{ii+1,i}, Mahald_data_cell{ii+9,i}] = fxn_nmt_Mahalanobis...
    (MPPCA_data_cell{ii, 5}.thrcov_PCA_score(Mahald_data_cell{2,reference_session_num},:), ... % reference
     MPPCA_data_cell{ii, 5}.thrcov_PCA_score(Mahald_data_cell{2,i},:));    % target
 
 % % Angle cal using all ca 1 bin data (normal calculation)
 [Angle_data_cell{ii+1,i}, Angle_data_cell{ii+9,i}] = fxn_nmt_cosine...
    (MPPCA_data_cell{ii, 3}(Mahald_data_cell{2,reference_session_num},:), ... % reference
     MPPCA_data_cell{ii, 3}(Mahald_data_cell{2,i},:));    % target
 
end
end

%% compare CS and ITI

CS_cells_960bin = MPPCA_data_cell{6, 3};
CS_cells_CS = mean(CS_cells_960bin(time_stamp_info{2,1},:),'all');
CS_cells_ITI1 = mean(CS_cells_960bin(time_stamp_info{4,1},:),'all');

CS_cells_data = {};
CS_cells_data_cell{1,1} = [CS_cells_CS,CS_cells_ITI1];

%% output
res_MPPCA_data_cell  = MPPCA_data_cell;
res_Angle_data_cell  = Angle_data_cell;
res_Mahald_data_cell = Mahald_data_cell;
res_CS_cells_data_cell = CS_cells_data_cell;

disp('MP-PCA-Mahalanobis cal finished!');
close all
%%
end
%%