%% CA1 connectivity analysis (210525 revised)
clc; clear; close all; disp('Previous data are cleared!'); 
%% determine connectivity calculating rate
% upsample_rate = 1; % for 1s binning data
% upsample_rate = 4; % for 200ms binning data, 4 x 50ms = 200ms.

% for 1s binning data
ca_binning_num = 20; % 20Hz to 1Hz data
timestamp_upsample_rate = 1; % upsaple with factor 1 (no upsample)

%% load mat-file
load('NatComm_data_CA1_calcium_data_raw_ctrl6_ko11.mat');

% ### select one variable. 
input_data = ca1_ctrl_m01;
% input_data = ca1_ctrl_m02;
% input_data = ca1_ctrl_m03;
% input_data = ca1_ctrl_m04;
% input_data = ca1_ctrl_m05;
% input_data = ca1_ctrl_m06;

% input_data = ca1_ko_m01;
% input_data = ca1_ko_m02;
% input_data = ca1_ko_m03;
% input_data = ca1_ko_m04;
% input_data = ca1_ko_m05;
% input_data = ca1_ko_m06;
% input_data = ca1_ko_m07;
% input_data = ca1_ko_m08;
% input_data = ca1_ko_m09;
% input_data = ca1_ko_m10;
% input_data = ca1_ko_m11;

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

end

%% output number of cells 
connectivity_results_num_cell = {};
connectivity_results_num_cell{1,1} = 'CS'; 
connectivity_results_num_cell{1,2} = 'US'; 
connectivity_results_num_cell{1,3} = 'LTM-CS'; 
connectivity_results_num_cell{1,4} = 'CS-US'; 
connectivity_results_num_cell{1,5} = 'CS-LTMCS'; 
connectivity_results_num_cell{1,6} = 'US-LTM'; 
connectivity_results_num_cell{1,7} = 'Triple'; 
connectivity_results_num_cell{1,8} = 'Total'; 

connectivity_results_num_cell{2,1} = numel(Cell_type_A_ID); % CS
connectivity_results_num_cell{2,2} = numel(Cell_type_B_ID); % US
connectivity_results_num_cell{2,3} = numel(Cell_type_C_ID); % Test-CS
connectivity_results_num_cell{2,4} = Cell_merge_A_B; % CS-US
connectivity_results_num_cell{2,5} = Cell_merge_A_C; % CS-TestCS
connectivity_results_num_cell{2,6} = Cell_merge_B_C; % US-Test-CS
connectivity_results_num_cell{2,7} = Cell_merbe_triple; % triple 
connectivity_results_num_cell{2,8} = size(G1_Ca_binned,2); % total

%% ##### Result output for fig 2c, 2d, and 2f. #####
result.num_cell      = connectivity_results_num_cell; % for fig 2c, 2d, 2f

%%