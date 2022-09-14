%% CA3 connectivity analysis (210528 revised)
function [result]  = fxn_CA3_triple_connectivity_cal_summary(input_data)
%%
% clc; clear; close all; disp('Previous data are cleared!'); tic;
%% determine connectivity calculating rate
% upsample_rate = 1; % for 1s binning data
% upsample_rate = 4; % for 500ms binning data, 4 x 50ms = 500ms.

% % for 1s binning data
% ca_binning_num = 20; % 20Hz to 1Hz data
% timestamp_upsample_rate = 1; % upsaple with factor 1 (no upsample)

% for 1s binning data
ca_binning_num = 10; % 20Hz to 10Hz data
timestamp_upsample_rate = 1; % upsaple with factor 1 (no upsample)

% for 500ms data
% ca_binning_num = 4; % 20Hz to 1Hz data
% timestamp_upsample_rate = 4;

%% filtering ans sigma-cut off
% ca_filt_data = fxn_filt_sigma_ca(ca_raw_data, sigma_if_z_ON);

sigma_cut_if_ON = 3; % dafault, Z score cutoff 
% sigma_cut_if_ON = 2; % Z score cutoff 
G1_concatenated_temp  = fxn_ca3_filt_sigma_cutoff(input_data, sigma_cut_if_ON);
%% Temporal binning
% [Ca_data_binned,  Ca_data_binned_mean] = funcHF_temporal_binning(Ca_data, binning_num)
binning_num = 20; % 20Hz to 1Hz data
[G1_Ca_binned,  G1_Ca_binned_mean] = fxn_temporal_binning(G1_concatenated_temp, binning_num);

%% input time stamp information
% upsample_rate = 1; % for 1s binning data
% upsample_rate = 4; % for 500ms binning data, 4 x 50ms = 500ms.

% normal cal.
time_stamp_info = fxn_ca3_ca55200_9sessions(timestamp_upsample_rate); % load session time information
disp(time_stamp_info)

% latter 9s ITI-1 ca.
% time_stamp_info = fxn_ca3_ca55200_9sessions_for_narrow_iti1(timestamp_upsample_rate); % load session time information

total_session_num  = 9;
responsiveness_threshold = 2;

G1_results = time_stamp_info;
%% Responsiveness calculation and ID-extracting
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
    
    % old code without isfinite
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
Cell_type_C_ID  = G1_results{8,5};  % {8 is LTM-CS+ cell-IDs for refference, fixed as 5}

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
    
    Cell_type_ABC_shared_ID   = intersect(Cell_type_A_B_shared_ID,Cell_type_C_ID);
    Cell_type_ABC_pure_A_ID   = setdiff(Cell_type_A_C_not_shared_A_ID, Cell_type_A_B_shared_ID);
    Cell_type_ABC_pure_B_ID   = setdiff(Cell_type_A_B_not_shared_B_ID, Cell_type_B_C_shared_ID);
    Cell_type_ABC_pure_C_ID   = setdiff(Cell_type_A_C_not_shared_C_ID, Cell_type_B_C_shared_ID);
%% count number of shared-cell
Cell_merge_A_B = numel(intersect(Cell_type_A_ID,Cell_type_B_ID));
Cell_merge_A_C = numel(intersect(Cell_type_A_ID,Cell_type_C_ID));
Cell_merge_B_C = numel(intersect(Cell_type_B_ID,Cell_type_C_ID));
Cell_merbe_triple = numel(intersect(intersect(Cell_type_A_ID,Cell_type_B_ID),Cell_type_C_ID));

Cell_merge_ITI1_A = numel(intersect(Cell_type_A_ID,Cell_type_ITI1_ID));
Cell_merge_ITI1_C = numel(intersect(Cell_type_ITI1_ID,Cell_type_C_ID));
Cell_merbe_ITI_triple = numel(intersect(intersect(Cell_type_A_ID, Cell_type_ITI1_ID),Cell_type_C_ID));

    Cell_type_ABC_shared_ID_num       = numel(Cell_type_ABC_shared_ID);
    Cell_type_ABC_pure_A_ID_num       = numel(Cell_type_ABC_pure_A_ID);
    Cell_type_ABC_pure_B_ID_num       = numel(Cell_type_ABC_pure_B_ID);
    Cell_type_ABC_pure_C_ID_num       = numel(Cell_type_ABC_pure_C_ID);
%% Calculation of Double occurence
% [raster_connectivity_t_index, raster_connectivity_index, i, ii] = ...
% fxn_sync_double_connectivity(data_input_temp, Cell_G1_ID, Cell_G2_ID, mymap_gradation_set, CS_ID, US_ID, LTMCS_ID)
% mymap_gradation_set, 0->Black, 1-> blue and green, 2-> red and green

% % A vs C., CS vs. LTMCS
% [raster_connectivity_t_index1, raster_connectivity_index1, i1, ii1] = ...
%     fxn_sync_double_connectivity(data_input_temp, Cell_type_A_ID, Cell_type_C_ID, 1, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('CS-LTMCS connectivity analysis done!');
% 
% % B vs C., US vs. LTMCS
% [raster_connectivity_t_index2, raster_connectivity_index2, i2, ii2] = ...
%     fxn_sync_double_connectivity(data_input_temp, Cell_type_B_ID, Cell_type_C_ID, 2, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('US-LTMCS connectivity analysis done!');
% 
% % A vs B., CS vs. US
% [raster_connectivity_t_index3, raster_connectivity_index3, i3, ii3] = ...
%     fxn_sync_double_connectivity(data_input_temp, Cell_type_A_ID, Cell_type_B_ID, 0, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('CS-LTMCS connectivity analysis done!');
% 
% display('Double connectivity calculation done!')
%% Calculation of Double occurence
% [raster_connectivity_t_index, raster_connectivity_index, i, ii] = ...
% fxn_sync_double_connectivity(data_input_temp, Cell_G1_ID, Cell_G2_ID, mymap_gradation_set, CS_ID, US_ID, LTMCS_ID)
% mymap_gradation_set, 0->Black, 1-> blue and green, 2-> red and green

% % A vs C., CS vs. LTMCS
% % Shared between CS and LTMCS
% [raster_connectivity_t_index_shared_AC, raster_connectivity_index_shared_AC, i1, ii1] = ...
%     fxn_sync_double_connectivity(data_input_temp, Cell_type_A_C_shared_ID, Cell_type_A_C_shared_ID, 1, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('CS-LTMCS mix done!');
% % Not-shared between CS and LTMCS (without shared subpopulation)
% [raster_connectivity_t_index_not_shared_AC, raster_connectivity_index_not_shared_AC, i, ii] = ...
%     fxn_sync_double_connectivity(data_input_temp, Cell_type_A_C_not_shared_A_ID, Cell_type_A_C_not_shared_C_ID, 2, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('CS-LTMCS shared done!');
% 
% % B vs C., US vs. LTMCS
% % Shared between CS and LTMCS
% [raster_connectivity_t_index_shared_BC, raster_connectivity_index_shared_BC, i1, ii1] = ...
%     fxn_sync_double_connectivity(data_input_temp, Cell_type_B_C_shared_ID, Cell_type_B_C_shared_ID, 1, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('US-LTMCS mix done!');
% % Not-shared between CS and LTMCS (without shared subpopulation)
% [raster_connectivity_t_index_not_shared_BC, raster_connectivity_index_not_shared_BC, i, ii] = ...
%     fxn_sync_double_connectivity(data_input_temp, Cell_type_B_C_not_shared_B_ID, Cell_type_B_C_not_shared_C_ID, 2, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('US-LTMCS shared done!');
% 
% % A vs B., CS vs. US
% % Shared between CS and LTMCS
% [raster_connectivity_t_index_shared_AB, raster_connectivity_index_shared_AB, i1, ii1] = ...
%     fxn_sync_double_connectivity(data_input_temp, Cell_type_A_B_shared_ID, Cell_type_A_B_shared_ID, 1, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('CS-US mix done!');
% % Not-shared between CS and LTMCS (without shared subpopulation)
% [raster_connectivity_t_index_not_shared_AB, raster_connectivity_index_not_shared_AB, i, ii] = ...
%     fxn_sync_double_connectivity(data_input_temp, Cell_type_A_B_not_shared_A_ID, Cell_type_A_B_not_shared_B_ID, 2, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('CS-US shared done!');
%% triple calculation

% triple
% [raster_connectivity_t_index_shared_triple, raster_connectivity_index_shared_triple, dim_Ai, dim_Aii, dim_Aiii] = ...
% fxn_sync_triple_connectivity(data_input_temp, Cell_type_ABC_shared_ID, Cell_type_ABC_shared_ID, Cell_type_ABC_shared_ID);
% disp('Triple shared done!');
% 
% [raster_connectivity_t_index_not_shared_triple, raster_connectivity_index_not_shared_triple, dim_Ai, dim_Aii, dim_Aiii] = ...
% fxn_sync_triple_connectivity(data_input_temp, Cell_type_ABC_pure_A_ID, Cell_type_ABC_pure_B_ID, Cell_type_ABC_pure_C_ID);
% disp('Triple not-shared done!');
% 
% display('Shared and not-shared subpopulation double connectivity calculation done!')
%% connectivity_summary
% % connectivity_table =[];
% % A vs C, CS vs. LTMCS % Shared between CS and LTMCS
% connectivity_results(1,i_index_cal) = raster_connectivity_index1;
% connectivity_results(2,i_index_cal) = raster_connectivity_index_shared_AC /2; % /2 -> because combination is calculated with same combination
% connectivity_results(3,i_index_cal) = raster_connectivity_index_not_shared_AC;
% % B vs C, US vs. LTMCS, % Shared between CS and LTMCS
% connectivity_results(4,i_index_cal) = raster_connectivity_index2;
% connectivity_results(5,i_index_cal) = raster_connectivity_index_shared_BC /2; % /2 -> because combination is calculated with same combination
% connectivity_results(6,i_index_cal) = raster_connectivity_index_not_shared_BC;
% % A vs B, CS vs. US, % Shared between CS and US
% connectivity_results(7,i_index_cal) = raster_connectivity_index3;
% connectivity_results(8,i_index_cal) = raster_connectivity_index_shared_AB /2; % /2 -> because combination is calculated with same combination
% connectivity_results(9,i_index_cal) = raster_connectivity_index_not_shared_AB;

% triple
% connectivity_results(10,i_index_cal) =  raster_connectivity_index_shared_triple /3;
% connectivity_results(11,i_index_cal) =  raster_connectivity_index_not_shared_triple;

%% connectivity_t_summary
% connectivity_results_t{1,i_index_cal+1} = raster_connectivity_t_index1;
% connectivity_results_t{2,i_index_cal+1} = raster_connectivity_t_index_shared_AC;
% connectivity_results_t{3,i_index_cal+1} = raster_connectivity_t_index_not_shared_AC;
% 
% connectivity_results_t{4,i_index_cal+1} = raster_connectivity_t_index2;
% connectivity_results_t{5,i_index_cal+1} = raster_connectivity_t_index_shared_BC;
% connectivity_results_t{6,i_index_cal+1} = raster_connectivity_t_index_not_shared_BC;
% 
% connectivity_results_t{7,i_index_cal+1} = raster_connectivity_t_index3;
% connectivity_results_t{8,i_index_cal+1} = raster_connectivity_t_index_shared_AB;
% connectivity_results_t{9,i_index_cal+1} = raster_connectivity_t_index_not_shared_AB;

% triple
% connectivity_results_t{10,i_index_cal+1} =  raster_connectivity_t_index_shared_triple /3;
% connectivity_results_t{11,i_index_cal+1} =  raster_connectivity_t_index_not_shared_triple;
%%
% loop_num = num2str(i_index_cal);
% disp_word_for_cal = ['Loop cal is now ', loop_num];
% disp(disp_word_for_cal);
%%
% close all; % to block figure 
% disp('Calulation finished! Use connectivity table for data analysis!')
end
%% for ITI 150s mean
%connectivity_results_t{x,5} = ('ITI-1'); connectivity_results_t{x,6} = ('ITI-2');
% for i_res = 1:9
% reshape_temp1 = reshape(connectivity_results_t{i_res,5}, 15, []);
% reshape_temp2 = reshape(connectivity_results_t{i_res,6}, 15, []);
% connectivity_results_t_block_ITI1(i_res,:) = mean(reshape_temp1 ,1);
% connectivity_results_t_block_ITI2(i_res,:) = mean(reshape_temp2 ,1);
% end
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

connectivity_results_num_cell{2,1} = numel(Cell_type_A_ID);
connectivity_results_num_cell{2,2} = numel(Cell_type_B_ID);
connectivity_results_num_cell{2,3} = numel(Cell_type_C_ID);
connectivity_results_num_cell{2,4} = Cell_merge_A_B;
connectivity_results_num_cell{2,5} = Cell_merge_A_C;
connectivity_results_num_cell{2,6} = Cell_merge_B_C;
connectivity_results_num_cell{2,7} = Cell_merbe_triple;
connectivity_results_num_cell{2,8} = size(G1_Ca_binned,2);

%% ########### 500ms code ###############
% ca_binning_num_500ms = 4; % 20Hz to 1Hz data
ca_binning_num_500ms = 10; % 500ms binning, 20Hz to 10Hz data

[G1_Ca_binned_500ms,  G1_Ca_binned_mean_500ms] = fxn_temporal_binning(G1_concatenated_temp, ca_binning_num_500ms);

% timestamp_upsample_rate_500ms = 4; % 50ms x 4 = 200
timestamp_upsample_rate_500ms = 2; % 1000ms / 500ms x  = 2
time_stamp_info_500ms = fxn_ca62400_12sessions(timestamp_upsample_rate_500ms); % load session time information

G1_results_500ms = time_stamp_info_500ms;
%% Responsiveness calculation and ID-extracting
for i = 1:total_session_num
G1_results_500ms{i,2} = G1_Ca_binned_500ms(time_stamp_info_500ms{i,1},:);
G1_results_500ms{i,3} = mean(G1_results_500ms{i,2},1);
G1_results_500ms{i,6} = G1_results_500ms{i,2} > 0;
end
%% Calculate responsiveness and cell-ID extraction
% find(X>0 & X<10,3)
for i = 2:total_session_num
G1_results_500ms{i,4} = (G1_results_500ms{i,3})./G1_results_500ms{1,3};
G1_results_500ms{i,5} = find(G1_results_500ms{i,4} >= responsiveness_threshold & isfinite(G1_results_500ms{i,4}));
G1_results_500ms{i,7} = find(isfinite(G1_results_500ms{i,4}));
    
    % old code without isfinite
    % G1_results_500ms{i,4} = (G1_results_500ms{i,3})./G1_results_500ms{1,3};
    % G1_results_500ms{i,5} = find(G1_results_500ms{i,4} >= responsiveness_threshold); % neuron ID
end
%% Synchronous activity 
% input parameter
connectivity_results_500ms_t = {};

connectivity_results_500ms_t{1,1} = ('CS vs LTM-CS'); 
connectivity_results_500ms_t{2,1} = ('US vs LTM-CS');
connectivity_results_500ms_t{3,1} = ('CS vs US');
connectivity_results_500ms_t{4,5} = ('ITI-1'); connectivity_results_500ms_t{4,6} = ('ITI-2');

connectivity_results_500ms_t{10,1} = ('Triple shared'); 
connectivity_results_500ms_t{11,1} = ('Triple not-shared');

for i_index_cal = 1:total_session_num
data_input_temp_500ms = G1_results_500ms{i_index_cal,6};  % {2 is input_data_session (2=CS, 3=US, 4=ITI), fixed as 6}

% Cell_type_A_ID  = G1_results_500ms{2,5};  % {2 is CS+ cell-IDs for refference, fixed as 5}
% Cell_type_B_ID  = G1_results_500ms{3,5};  % {3 is US+ cell-IDs for refference, fixed as 5}
% Cell_type_ITI1_ID  = G1_results_500ms{4,5};  % {3 is US+ cell-IDs for refference, fixed as 5}
% Cell_type_ITI2_ID  = G1_results_500ms{5,5};  % {3 is US+ cell-IDs for refference, fixed as 5}
% Cell_type_C_ID  = G1_results_500ms{8,5};  % {8 is LTM-CS+ cell-IDs for refference, fixed as 5}
% 
%     % 210524 added
%     Cell_type_A_C_shared_ID   = intersect(Cell_type_A_ID,Cell_type_C_ID);
%     Cell_type_A_C_not_shared_A_ID = setdiff(Cell_type_A_ID, Cell_type_A_C_shared_ID);
%     Cell_type_A_C_not_shared_C_ID = setdiff(Cell_type_C_ID, Cell_type_A_C_shared_ID);
% 
%     Cell_type_A_B_shared_ID   = intersect(Cell_type_A_ID,Cell_type_B_ID);
%     Cell_type_A_B_not_shared_A_ID = setdiff(Cell_type_A_ID, Cell_type_A_B_shared_ID);
%     Cell_type_A_B_not_shared_B_ID = setdiff(Cell_type_B_ID, Cell_type_A_B_shared_ID);
% 
%     Cell_type_B_C_shared_ID   = intersect(Cell_type_B_ID,Cell_type_C_ID);
%     Cell_type_B_C_not_shared_B_ID = setdiff(Cell_type_B_ID, Cell_type_B_C_shared_ID);
%     Cell_type_B_C_not_shared_C_ID = setdiff(Cell_type_C_ID, Cell_type_B_C_shared_ID);
%     
%     Cell_type_A_B_C_shared_ID   = intersect(Cell_type_A_B_shared_ID,Cell_type_C_ID);
% 
%% count number of shared-cell
% Cell_merge_A_B = numel(intersect(Cell_type_A_ID,Cell_type_B_ID));
% Cell_merge_A_C = numel(intersect(Cell_type_A_ID,Cell_type_C_ID));
% Cell_merge_B_C = numel(intersect(Cell_type_B_ID,Cell_type_C_ID));
% Cell_merbe_triple = numel(intersect(intersect(Cell_type_A_ID,Cell_type_B_ID),Cell_type_C_ID));
% 
% Cell_merge_ITI1_A = numel(intersect(Cell_type_A_ID,Cell_type_ITI1_ID));
% Cell_merge_ITI1_C = numel(intersect(Cell_type_ITI1_ID,Cell_type_C_ID));
% Cell_merbe_ITI_triple = numel(intersect(intersect(Cell_type_A_ID, Cell_type_ITI1_ID),Cell_type_C_ID));
%% Calculation of Double occurence
% [raster_connectivity_t_index, raster_connectivity_index, i, ii] = ...
% fxn_sync_double_connectivity(data_input_temp_500ms, Cell_G1_ID, Cell_G2_ID, mymap_gradation_set, CS_ID, US_ID, LTMCS_ID)
% mymap_gradation_set, 0->Black, 1-> blue and green, 2-> red and green

% % A vs C., CS vs. LTMCS
% [raster_connectivity_500ms_t_index1, raster_connectivity_500ms_index1, i1, ii1] = ...
%     fxn_sync_double_connectivity(data_input_temp_500ms, Cell_type_A_ID, Cell_type_C_ID, 1, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('CS-LTMCS connectivity analysis done!');
% 
% % B vs C., US vs. LTMCS
% [raster_connectivity_500ms_t_index2, raster_connectivity_500ms_index2, i2, ii2] = ...
%     fxn_sync_double_connectivity(data_input_temp_500ms, Cell_type_B_ID, Cell_type_C_ID, 2, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('US-LTMCS connectivity analysis done!');
% 
% % A vs B., CS vs. US
% [raster_connectivity_500ms_t_index3, raster_connectivity_500ms_index3, i3, ii3] = ...
%     fxn_sync_double_connectivity(data_input_temp_500ms, Cell_type_A_ID, Cell_type_B_ID, 0, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('CS-LTMCS connectivity analysis done!');
% 
% display('Double connectivity calculation done!')
%% Calculation of Double occurence
% [raster_connectivity_500ms_t_index, raster_connectivity_500ms_index, i, ii] = ...
% fxn_sync_double_connectivity(data_input_temp_500ms, Cell_G1_ID, Cell_G2_ID, mymap_gradation_set, CS_ID, US_ID, LTMCS_ID)
% mymap_gradation_set, 0->Black, 1-> blue and green, 2-> red and green

% % A vs C., CS vs. LTMCS
% % Shared between CS and LTMCS
% [raster_connectivity_500ms_t_index_shared_AC, raster_connectivity_500ms_index_shared_AC, i1, ii1] = ...
%     fxn_sync_double_connectivity(data_input_temp_500ms, Cell_type_A_C_shared_ID, Cell_type_A_C_shared_ID, 1, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('CS-LTMCS mix done!');
% % Not-shared between CS and LTMCS (without shared subpopulation)
% [raster_connectivity_500ms_t_index_not_shared_AC, raster_connectivity_500ms_index_not_shared_AC, i, ii] = ...
%     fxn_sync_double_connectivity(data_input_temp_500ms, Cell_type_A_C_not_shared_A_ID, Cell_type_A_C_not_shared_C_ID, 2, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('CS-LTMCS shared done!');
% 
% % B vs C., US vs. LTMCS
% % Shared between CS and LTMCS
% [raster_connectivity_500ms_t_index_shared_BC, raster_connectivity_500ms_index_shared_BC, i1, ii1] = ...
%     fxn_sync_double_connectivity(data_input_temp_500ms, Cell_type_B_C_shared_ID, Cell_type_B_C_shared_ID, 1, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('US-LTMCS mix done!');
% % Not-shared between CS and LTMCS (without shared subpopulation)
% [raster_connectivity_500ms_t_index_not_shared_BC, raster_connectivity_500ms_index_not_shared_BC, i, ii] = ...
%     fxn_sync_double_connectivity(data_input_temp_500ms, Cell_type_B_C_not_shared_B_ID, Cell_type_B_C_not_shared_C_ID, 2, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('US-LTMCS shared done!');
% 
% % A vs B., CS vs. US
% % Shared between CS and LTMCS
% [raster_connectivity_500ms_t_index_shared_AB, raster_connectivity_500ms_index_shared_AB, i1, ii1] = ...
%     fxn_sync_double_connectivity(data_input_temp_500ms, Cell_type_A_B_shared_ID, Cell_type_A_B_shared_ID, 1, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('CS-US mix done!');
% % Not-shared between CS and LTMCS (without shared subpopulation)
% [raster_connectivity_500ms_t_index_not_shared_AB, raster_connectivity_500ms_index_not_shared_AB, i, ii] = ...
%     fxn_sync_double_connectivity(data_input_temp_500ms, Cell_type_A_B_not_shared_A_ID, Cell_type_A_B_not_shared_B_ID, 2, ...
%     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID); disp('CS-US shared done!');
% 
% display('Double shared and not-shared subpopulation double connectivity calculation done!')
%% triple calculation

% triple
[raster_connectivity_500ms_t_index_shared_triple, raster_connectivity_500ms_index_shared_triple, dim_Ai, dim_Aii, dim_Aiii,shared_total_sync, shared_cell_500ms, shared_time_500ms, ...
    shared_raster_triple_connectivity_raw, shared_standard_connectivity] = ...
fxn_sync_triple_connectivity(data_input_temp_500ms, Cell_type_ABC_shared_ID, Cell_type_ABC_shared_ID, Cell_type_ABC_shared_ID);
disp('Triple shared done!');

[raster_connectivity_500ms_t_index_not_shared_triple, raster_connectivity_500ms_index_not_shared_triple, dim_Ai, dim_Aii, dim_Aiii,non_shared_total_sync, non_shared_cell_500ms, non_shared_time_500ms, ...
    non_shared_raster_triple_connectivity_raw, non_shared_standard_connectivity] = ...
fxn_sync_triple_connectivity(data_input_temp_500ms, Cell_type_ABC_pure_A_ID, Cell_type_ABC_pure_B_ID, Cell_type_ABC_pure_C_ID);
disp('Triple not-shared done!');

display('Triple shared and not-shared subpopulation double connectivity calculation done!')
%% connectivity_summary
% connectivity_table =[];
% % A vs C, CS vs. LTMCS % Shared between CS and LTMCS
% connectivity_results_500ms(1,i_index_cal) = raster_connectivity_500ms_index1;
% connectivity_results_500ms(2,i_index_cal) = raster_connectivity_500ms_index_shared_AC /2; % /2 -> because combination is calculated with same combination
% connectivity_results_500ms(3,i_index_cal) = raster_connectivity_500ms_index_not_shared_AC;
% % B vs C, US vs. LTMCS, % Shared between CS and LTMCS
% connectivity_results_500ms(4,i_index_cal) = raster_connectivity_500ms_index2;
% connectivity_results_500ms(5,i_index_cal) = raster_connectivity_500ms_index_shared_BC /2; % /2 -> because combination is calculated with same combination
% connectivity_results_500ms(6,i_index_cal) = raster_connectivity_500ms_index_not_shared_BC;
% % A vs B, CS vs. US, % Shared between CS and US
% connectivity_results_500ms(7,i_index_cal) = raster_connectivity_500ms_index3;
% connectivity_results_500ms(8,i_index_cal) = raster_connectivity_500ms_index_shared_AB /2; % /2 -> because combination is calculated with same combination
% connectivity_results_500ms(9,i_index_cal) = raster_connectivity_500ms_index_not_shared_AB;

% triple
connectivity_results_500ms(10,i_index_cal) =  raster_connectivity_500ms_index_shared_triple /3;
connectivity_results_500ms(11,i_index_cal) =  raster_connectivity_500ms_index_not_shared_triple;

result.total_sync(10,i_index_cal) = shared_total_sync;
result.total_sync(11,i_index_cal) = non_shared_total_sync;

result.total_cell_500ms(10,i_index_cal) = shared_cell_500ms;
result.total_cell_500ms(11,i_index_cal) = non_shared_cell_500ms;

result.time_500ms(10,i_index_cal) = shared_time_500ms;
result.time_500ms(11,i_index_cal) = non_shared_time_500ms;

result.triple_sync_raw_500ms(10,i_index_cal) = shared_raster_triple_connectivity_raw;
result.triple_sync_raw_500ms(11,i_index_cal) = non_shared_raster_triple_connectivity_raw;

result.triple_sync_standard_500ms(10,i_index_cal) = shared_standard_connectivity;
result.triple_sync_standard_500ms(11,i_index_cal) = non_shared_standard_connectivity;

%% connectivity_t_summary
% connectivity_results_500ms_t{1,i_index_cal+1} = raster_connectivity_500ms_t_index1;
% connectivity_results_500ms_t{2,i_index_cal+1} = raster_connectivity_500ms_t_index_shared_AC;
% connectivity_results_500ms_t{3,i_index_cal+1} = raster_connectivity_500ms_t_index_not_shared_AC;
% 
% connectivity_results_500ms_t{4,i_index_cal+1} = raster_connectivity_500ms_t_index2;
% connectivity_results_500ms_t{5,i_index_cal+1} = raster_connectivity_500ms_t_index_shared_BC;
% connectivity_results_500ms_t{6,i_index_cal+1} = raster_connectivity_500ms_t_index_not_shared_BC;
% 
% connectivity_results_500ms_t{7,i_index_cal+1} = raster_connectivity_500ms_t_index3;
% connectivity_results_500ms_t{8,i_index_cal+1} = raster_connectivity_500ms_t_index_shared_AB;
% connectivity_results_500ms_t{9,i_index_cal+1} = raster_connectivity_500ms_t_index_not_shared_AB;

% triple
connectivity_results_500ms_t{10,i_index_cal+1} =  raster_connectivity_500ms_t_index_shared_triple /3;
connectivity_results_500ms_t{11,i_index_cal+1} =  raster_connectivity_500ms_t_index_not_shared_triple;
%%
loop_num = num2str(i_index_cal);
disp_word_for_cal = ['Loop cal for 500ms is now ', loop_num];
disp(disp_word_for_cal);
%%
close all; % to block figure 
disp('Calulation for 500ms version finished! Use connectivity table for data analysis!')
end
%% for ITI 150s mean
%connectivity_results_500ms_t{x,5} = ('ITI-1'); connectivity_results_500ms_t{x,6} = ('ITI-2');
% for i_res = 1:11 % change from 1:9 to 1:11
% reshape_temp1_500ms = reshape(connectivity_results_500ms_t{i_res,5}, 15*timestamp_upsample_rate_500ms, []);
% reshape_temp2_500ms = reshape(connectivity_results_500ms_t{i_res,6}, 15*timestamp_upsample_rate_500ms, []);
% connectivity_results_500ms_t_block_ITI1(i_res,:) = mean(reshape_temp1_500ms ,1);
% connectivity_results_500ms_t_block_ITI2(i_res,:) = mean(reshape_temp2_500ms ,1);
% end
%% output number of cells 
connectivity_results_500ms_num_cell = {};
connectivity_results_500ms_num_cell{1,1} = 'CS'; 
connectivity_results_500ms_num_cell{1,2} = 'US'; 
connectivity_results_500ms_num_cell{1,3} = 'LTM-CS'; 
connectivity_results_500ms_num_cell{1,4} = 'CS-US'; 
connectivity_results_500ms_num_cell{1,5} = 'CS-LTMCS'; 
connectivity_results_500ms_num_cell{1,6} = 'US-LTM'; 
connectivity_results_500ms_num_cell{1,7} = 'Triple'; 
connectivity_results_500ms_num_cell{1,8} = 'Total'; 

connectivity_results_500ms_num_cell{2,1} = numel(Cell_type_A_ID);
connectivity_results_500ms_num_cell{2,2} = numel(Cell_type_B_ID);
connectivity_results_500ms_num_cell{2,3} = numel(Cell_type_C_ID);
connectivity_results_500ms_num_cell{2,4} = Cell_merge_A_B;
connectivity_results_500ms_num_cell{2,5} = Cell_merge_A_C;
connectivity_results_500ms_num_cell{2,6} = Cell_merge_B_C;
connectivity_results_500ms_num_cell{2,7} = Cell_merbe_triple;
connectivity_results_500ms_num_cell{2,8} = size(G1_Ca_binned,2);
%% ########## 500 ms analysis part finished #########
disp('500 ms analysis part finished.')
%% figure connectivity results_t 
figure; 
% subplot(221); imagesc(connectivity_results);
% subplot(223); imagesc(connectivity_results_t_block_ITI1);
subplot(222); imagesc(connectivity_results_500ms);
% subplot(224); imagesc(connectivity_results_500ms_t_block_ITI1);

%% Result output
result.num_cell      = connectivity_results_num_cell;

% result.summary       = connectivity_results;
% result.ITI1          = connectivity_results_t_block_ITI1;
% result.ITI2          = connectivity_results_t_block_ITI2;

result.summary_500ms = connectivity_results_500ms;
% result.ITI1_500ms    = connectivity_results_500ms_t_block_ITI1;
% result.ITI2_500ms    = connectivity_results_500ms_t_block_ITI2;
%%

%%
end
