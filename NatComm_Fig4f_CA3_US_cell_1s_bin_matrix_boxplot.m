%% Responsiveness analysis
clc; clear; close all; disp('Previous data are cleared!')
%% load mat-file
load("NatComm_data_head_fixed_ca1_ca3_data_set.mat");

%% generate concatenated data

% #### CA1-HF data ####
% G1_concatenated = [hf_ca1_ctrl_m01, hf_ca1_ctrl_m02]; % for paper           
% G2_concatenated = [hf_ca1_ko_m01, hf_ca1_ko_m02]; % for paper

% #### CA3-HF data ####
G1_concatenated = [hf_ca3_ctrl_m01, hf_ca3_ctrl_m02, hf_ca3_ctrl_m03, hf_ca3_ctrl_m04]; % for paper     
G2_concatenated = [hf_ca3_ko_m01, hf_ca3_ko_m02]; % for paper     

%% filtering ans sigma-cut off
% ca_filt_data = fxn_filt_sigma_ca(ca_raw_data, sigma_if_z_ON);

sigma_cut_if_ON = 3; % Z score cutoff 

G1_concatenated_temp  = fxnHF_filt_sigma_cutoff(G1_concatenated, sigma_cut_if_ON);
G2_concatenated_temp  = fxnHF_filt_sigma_cutoff(G2_concatenated, sigma_cut_if_ON);
%% Temporal binning
% [Ca_data_binned,  Ca_data_binned_mean] = funcHF_temporal_binning(Ca_data, binning_num)

binning_num = 20; % 20Hz to 1Hz data

[G1_Ca_binned,  G1_Ca_binned_mean] = fxnHF_temporal_binning(G1_concatenated_temp, binning_num);
[G2_Ca_binned,  G2_Ca_binned_mean] = fxnHF_temporal_binning(G2_concatenated_temp, binning_num);

reshape_range1 = (376:960); % for CS
reshape_range2 = (976:1560); % for US
reshape_range3 = (1576:2160); % for extra

% fxnHF_1s_bin_matrix(G1_Ca_binned, G2_Ca_binned, G1_Ca_binned_mean, G2_Ca_binned_mean, reshape_range1)  
% fxnHF_1s_bin_matrix(G1_Ca_binned, G2_Ca_binned, G1_Ca_binned_mean, G2_Ca_binned_mean, reshape_range2)  
% fxnHF_1s_bin_matrix(G1_Ca_binned, G2_Ca_binned, G1_Ca_binned_mean, G2_Ca_binned_mean, reshape_range3)   
%% input time stamp information
% onset_bin1 = 360; % only CS session
onset_bin1 = 960; % only US session
% onset_bin1 = 1560; % additional session

results = fxnHF_ca57600_8sessions(onset_bin1); % load session time information
total_session_num  = 8;
responsiveness_threshold = 2;

G1_results = results;
G2_results = results;

%% Calculate mean in each session
for i = 1:total_session_num
G1_results{i,2} = G1_Ca_binned(results{i,1},:);
G1_results{i,3} = mean(G1_results{i,2},1);
G1_results{1,6} = G1_Ca_binned./G1_results{1,3};

G2_results{i,2} = G2_Ca_binned(results{i,1},:);
G2_results{i,3} = mean(G2_results{i,2},1);
G2_results{1,6} = G2_Ca_binned./G2_results{1,3};
end
%% Calculate responsiveness and cell-ID extraction
% find(X>0 & X<10,3)
for i = 2:total_session_num
G1_results{i,4} = (G1_results{i,3})./G1_results{1,3};
G1_results{i,5} = find(G1_results{i,4} >= responsiveness_threshold & isfinite(G1_results{i,4}));
G1_results{i,7} = find(isfinite(G1_results{i,4}));

G2_results{i,4} = (G2_results{i,3})./G2_results{1,3};
% G2_results{i,5} = find(G2_results{i,4} >= responsiveness_threshold); old code
G2_results{i,5} = find(G2_results{i,4} >= responsiveness_threshold & isfinite(G2_results{i,4}));
G2_results{i,7} = find(isfinite(G2_results{i,4}));
end
%% pick cell ID

Cell_type_A_ID_G1  = G1_results{2,5};  % {2 is CS+ cell-IDs for refference, fixed as 5}
Cell_type_B_ID_G1  = G1_results{3,5};  % {3 is US+ cell-IDs for refference, fixed as 5}
Cell_type_C_ID_G1  = G1_results{7,5}; % {13 is LTM-CS+ cell-IDs for refference, fixed as 5} % ltm-cs
Cell_type_ITI1_ID_G1  = G1_results{4,5};  % {3 is US+ cell-IDs for refference, fixed as 5}
Cell_type_ITI2_ID_G1  = G1_results{5,5};  % {3 is US+ cell-IDs for refference, fixed as 5}

Cell_type_A_ID_G2  = G2_results{2,5};  % {2 is CS+ cell-IDs for refference, fixed as 5}
Cell_type_B_ID_G2  = G2_results{3,5};  % {3 is US+ cell-IDs for refference, fixed as 5}
Cell_type_C_ID_G2  = G2_results{7,5}; % {13 is LTM-CS+ cell-IDs for refference, fixed as 5} % ltm-cs
Cell_type_ITI1_ID_G2  = G2_results{4,5};  % {3 is US+ cell-IDs for refference, fixed as 5}
Cell_type_ITI2_ID_G2  = G2_results{5,5};  % {3 is US+ cell-IDs for refference, fixed as 5}


%% Calculate reactivation in responsiveness-classified groups
for i = 2:total_session_num
    for ii = 1:total_session_num
        G1_results{15+ii,i} =  G1_results{ii,2}(:,G1_results{i,5});
        G1_results{30+ii,i} =  mean(G1_results{ii,2}(:,G1_results{i,5}),1);
        G1_results{45+ii,i} =  mean(G1_results{30+ii,i},2);
        
        G2_results{15+ii,i} =  G2_results{ii,2}(:,G2_results{i,5});
        G2_results{30+ii,i} =  mean(G2_results{ii,2}(:,G2_results{i,5}),1);
        G2_results{45+ii,i} =  mean(G2_results{30+ii,i},2);
    end
end
%% Calculate Responsiveness matrix (mean, error, and ranksum-stat results p h)
for i = 2:total_session_num
    for ii = 2:total_session_num
        G1_results_responsive{ii,i} =  G1_results{ii,4}(:,G1_results{i,5});
        G1_results_responsive{15+ii,i} =  mean(G1_results_responsive{ii,i},2);
        G1_results_responsive{30+ii,i} =  std(G1_results_responsive{ii,i})./sqrt(numel(G1_results_responsive{ii,i}));
        
        G2_results_responsive{ii,i} =  G2_results{ii,4}(:,G2_results{i,5});
        G2_results_responsive{15+ii,i} =  mean(G2_results_responsive{ii,i},2);
        G2_results_responsive{30+ii,i} =  std(G2_results_responsive{ii,i})./sqrt(numel(G1_results_responsive{ii,i}));
        
        [G1_G2_results_responsive_ranksum_p(ii-1,i-1), G1_G2_results_responsive_ranksum_h(ii-1,i-1)] = ranksum((G1_results_responsive{ii,i}),(G2_results_responsive{ii,i}));     
        
    end
end
        G1_results_responsive_mean  = cell2mat(G1_results_responsive(17:23,2:8)); % mean
        G1_results_responsive_error = cell2mat(G1_results_responsive(32:38,2:8)); % error
        
        G2_results_responsive_mean  = cell2mat(G2_results_responsive(17:23,2:8)); % mean
        G2_results_responsive_error = cell2mat(G2_results_responsive(32:38,2:8)); % error
        
        G1_G2_results_responsive_subtracttion = G1_results_responsive_mean-G2_results_responsive_mean;
%% mean, error, and ranksum-stat figure
% figure;
% subplot(232); imagesc(G1_results_responsive_mean)
% subplot(233); imagesc(G2_results_responsive_mean)
% subplot(234); imagesc(G1_G2_results_responsive_subtracttion)
% subplot(235); imagesc(G1_G2_results_responsive_ranksum_p)
% subplot(236); imagesc(G1_G2_results_responsive_ranksum_h)
%% Baseline-responsiveness zero-cell removed Reponsive 1s bin matrix (ALL)
% (ALL cells resutls, removed IDs showeing infinite value, due to baseline zero value)

% reshape_range   = (16:600);
% reshape_range = (151:600);

G1_Ca_binned_responsive = G1_results{1,6}(:,G1_results{2,7}); G1_Ca_binned_responsive_mean = mean(G1_Ca_binned_responsive,2);
G2_Ca_binned_responsive = G2_results{1,6}(:,G2_results{2,7}); G2_Ca_binned_responsive_mean = mean(G2_Ca_binned_responsive,2);

% fxnHF_1s_bin_matrix(G1_Ca_binned_responsive, G2_Ca_binned_responsive, ...
%                       G1_Ca_binned_responsive_mean, G2_Ca_binned_responsive_mean, reshape_range1)
% fxnHF_1s_bin_matrix(G1_Ca_binned_responsive, G2_Ca_binned_responsive, ...
%                       G1_Ca_binned_responsive_mean, G2_Ca_binned_responsive_mean, reshape_range2)
%% Baseline-responsiveness zero-cell removed Reponsive 1s bin matrix (CS,US, and LTM-CS)
% (CS,US, and LTM-CS cells resutls, removed IDs showeing infinite value, due to baseline zero value)
% CS, US, LTM-CS cell pickup

% reshape_range   = (16:600);
% reshape_range = (151:600);

% G1_neruon_ID = unique([G1_results{2,5}, G1_results{3,5}, G1_results{7,5}]);
% G2_neruon_ID = unique([G2_results{2,5}, G2_results{3,5}, G2_results{7,5}]);
% 
% G1_Ca_binned_responsive = G1_results{1,6}(:,G1_neruon_ID); G1_Ca_binned_responsive_mean = mean(G1_Ca_binned_responsive,2);
% G2_Ca_binned_responsive = G2_results{1,6}(:,G2_neruon_ID); G2_Ca_binned_responsive_mean = mean(G2_Ca_binned_responsive,2);
% 
% fxnHF_1s_bin_matrix(G1_Ca_binned_responsive, G2_Ca_binned_responsive, ...
%                       G1_Ca_binned_responsive_mean, G2_Ca_binned_responsive_mean, reshape_range1)   
% fxnHF_1s_bin_matrix(G1_Ca_binned_responsive, G2_Ca_binned_responsive, ...
%                       G1_Ca_binned_responsive_mean, G2_Ca_binned_responsive_mean, reshape_range2)   
%% Box plot
% fxn_boxplot_double(G1_results_ratio, G2_results_ratio, input_session_num, section_num);
% section_num_end = 8; % con-ko -> 14
% 
% for i = 2:section_num_end
%  fxnHF_boxplot_double(G1_results_responsive, G2_results_responsive, total_session_num, i);
% end
%%
close all
%% ##### Z-score data using ####


%% Baseline-responsiveness zero-cell removed Z-score 1s bin matrix
% (ALL cells resutls, removed IDs showeing infinite value, due to baseline zero value)

% reshape_range   = (16:600);
% reshape_range = (151:600);

% G1_Ca_binned_resp0_removed_Z = G1_Ca_binned(:,G1_results{2,7}); G1_Ca_binned_resp0_removed_Z_mean = mean(G1_Ca_binned_resp0_removed_Z,2);
% G2_Ca_binned_resp0_removed_Z = G2_Ca_binned(:,G2_results{2,7}); G2_Ca_binned_resp0_removed_Z_mean = mean(G2_Ca_binned_resp0_removed_Z,2);
% 
% fxnHF_1s_bin_matrix(G1_Ca_binned_resp0_removed_Z, G2_Ca_binned_resp0_removed_Z, ...
%                       G1_Ca_binned_resp0_removed_Z_mean, G2_Ca_binned_resp0_removed_Z_mean, reshape_range1)
% fxnHF_1s_bin_matrix(G1_Ca_binned_resp0_removed_Z, G2_Ca_binned_resp0_removed_Z, ...
%                       G1_Ca_binned_resp0_removed_Z_mean, G2_Ca_binned_resp0_removed_Z_mean, reshape_range2)
%% Baseline-responsiveness zero-cell removed Z-score matrix (mean, error, and ranksum-stat results p h)
for i = 2:total_session_num
    for ii = 1:total_session_num
        G1_results_resp0_removed_Z{ii,i} =  G1_results{ii,3}(:,G1_results{i,5});
        G1_results_resp0_removed_Z{15+ii,i} =  mean(G1_results_resp0_removed_Z{ii,i},2);
        G1_results_resp0_removed_Z{30+ii,i} =  std(G1_results_resp0_removed_Z{ii,i})./sqrt(numel(G1_results_resp0_removed_Z{ii,i}));
        
        G2_results_resp0_removed_Z{ii,i} =  G2_results{ii,3}(:,G2_results{i,5});
        G2_results_resp0_removed_Z{15+ii,i} =  mean(G2_results_resp0_removed_Z{ii,i},2);
        G2_results_resp0_removed_Z{30+ii,i} =  std(G2_results_resp0_removed_Z{ii,i})./sqrt(numel(G1_results_resp0_removed_Z{ii,i}));
        
        [G1_G2_results_resp0_removed_Z_ranksum_p(ii,i), G1_G2_results_resp0_removed_Z_ranksum_h(ii,i) ...
            G1_G2_results_resp0_removed_Z_ranksum_stats(ii,i)] = ranksum((G1_results_resp0_removed_Z{ii,i}),(G2_results_resp0_removed_Z{ii,i}));           
            G1_G2_results_resp0_removed_Z_ranksum_stats_zval(ii,i) = G1_G2_results_resp0_removed_Z_ranksum_stats(ii, i).zval;
        
    end
end
        G1_results_resp0_removed_Z_mean  = cell2mat(G1_results_resp0_removed_Z(16:23,2:8)); % mean
        G1_results_resp0_removed_Z_error = cell2mat(G1_results_resp0_removed_Z(31:38,2:8)); % error
        
        G2_results_resp0_removed_Z_mean  = cell2mat(G2_results_resp0_removed_Z(16:23,2:8)); % mean
        G2_results_resp0_removed_Z_error = cell2mat(G2_results_resp0_removed_Z(31:38,2:8)); % error
        
        G1_G2_results_resp0_removed_Z_subtracttion = G1_results_resp0_removed_Z_mean-G2_results_resp0_removed_Z_mean;
%% Baseline-responsiveness zero-cell removed Z-score 1s bin matrix  
% (CS,US, and LTM-CS cells resutls, removed IDs showeing infinite value, due to baseline zero value)
% CS, US, LTM-CS cell pickup

reshape_range   = (16:600);
% reshape_range = (151:600);

% G1_neruon_ID = unique([G1_results{2,5}]); % CS
% G2_neruon_ID = unique([G2_results{2,5}]); % CS

G1_neruon_ID = unique([G1_results{3,5}]); % US
G2_neruon_ID = unique([G2_results{3,5}]); % US

% G1_neruon_ID = unique([G1_results{7,5}]); % LTM-CS
% G2_neruon_ID = unique([G2_results{7,5}]); % LTM-CS

% G1_neruon_ID = unique([G1_results{2,5}, G1_results{3,5}, G1_results{7,5}]); % CS,US, and LTM-CS
% G2_neruon_ID = unique([G2_results{2,5}, G2_results{3,5}, G2_results{7,5}]); % CS,US, and LTM-CS

G1_Ca_binned_resp0_removed_Z = G1_Ca_binned(:,G1_neruon_ID); G1_Ca_binned_resp0_removed_Z_mean = mean(G1_Ca_binned_resp0_removed_Z,2);
G2_Ca_binned_resp0_removed_Z = G2_Ca_binned(:,G2_neruon_ID); G2_Ca_binned_resp0_removed_Z_mean = mean(G2_Ca_binned_resp0_removed_Z,2);

% fxnHF_1s_bin_matrix(G1_Ca_binned_resp0_removed_Z, G2_Ca_binned_resp0_removed_Z, ...
%                       G1_Ca_binned_resp0_removed_Z_mean, G2_Ca_binned_resp0_removed_Z_mean, reshape_range1)   
fxnHF_1s_bin_matrix(G1_Ca_binned_resp0_removed_Z, G2_Ca_binned_resp0_removed_Z, ...
                      G1_Ca_binned_resp0_removed_Z_mean, G2_Ca_binned_resp0_removed_Z_mean, reshape_range2)   
% fxnHF_1s_bin_matrix(G1_Ca_binned_resp0_removed_Z, G2_Ca_binned_resp0_removed_Z, ...
%                       G1_Ca_binned_resp0_removed_Z_mean, G2_Ca_binned_resp0_removed_Z_mean, reshape_range3)   

%% mean, error, and ranksum-stat figure
figure;
subplot(232); imagesc(G1_results_resp0_removed_Z_mean)
subplot(233); imagesc(G2_results_resp0_removed_Z_mean)
subplot(234); imagesc(G1_G2_results_resp0_removed_Z_subtracttion)
subplot(235); imagesc(G1_G2_results_resp0_removed_Z_ranksum_p)
subplot(236); imagesc(G1_G2_results_resp0_removed_Z_ranksum_h)
%% Box plot for CS-only
% fxn_boxplot_double(G1_results_ratio, G2_results_ratio, input_session_num, section_num);
% boxplot_section_num_end = 4; % 9=STM removed,  12=Rest summarised, 14=entire
% 
% G1_results_resp0_removed_Z_CS_only = G1_results_resp0_removed_Z([1:2 4:5],[1:2 4:5]);
% G2_results_resp0_removed_Z_CS_only = G2_results_resp0_removed_Z([1:2 4:5],[1:2 4:5]);
% box_name_cell = {'Baseline', 'CS','ITI-1','ITI-2'};
% 
% for i = 1:boxplot_section_num_end
%  fxnHF_boxplot_double(G1_results_resp0_removed_Z_CS_only, G2_results_resp0_removed_Z_CS_only, boxplot_section_num_end, i,box_name_cell);
% end
%% Box plot for US-only
% fxn_boxplot_double(G1_results_ratio, G2_results_ratio, input_session_num, section_num);
boxplot_section_num_end = 4; % 9=STM removed,  12=Rest summarised, 14=entire

G1_results_resp0_removed_Z_US_only = G1_results_resp0_removed_Z([1 3:5],[1 3:5]);
G2_results_resp0_removed_Z_US_only = G2_results_resp0_removed_Z([1 3:5],[1 3:5]);
box_name_cell = {'Baseline', 'US','ITI-1','ITI-2'};

for i = 1:boxplot_section_num_end
 fxnHF_boxplot_double(G1_results_resp0_removed_Z_US_only, G2_results_resp0_removed_Z_US_only, boxplot_section_num_end, i,box_name_cell);
end
%%
%% figure Color-raster plot for entire activity
% fxn_plot_CS_US_LTM_cells(data_input_temp, Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID)

% fxn_plot_CS_US_LTM_cells(logical(G1_Ca_binned), Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID)

% %% function code for figure
% temp_cell_id_G1 = Cell_type_A_ID_G1; % CS cell
% temp_cell_id = Cell_type_B_ID; % US cell
% temp_cell_id = Cell_type_C_ID; % LTM-CS cell

% temp_cell_id_G2 = Cell_type_A_ID_G2; % CS cell
% temp_cell_id = Cell_type_B_ID; % US cell
% temp_cell_id = Cell_type_C_ID; % LTM-CS cell

% reshape_range1 = (376:960); % for CS
% reshape_range2 = (976:1560); % for US
% reshape_range3 = (1576:2160); % for CSUS

cell_id_G1      = temp_cell_id_G1;
cell_id_G2      = temp_cell_id_G2;


temp_range1 = [361:960];
figure_range1 = [90:600];

% temp_range2 = [2760:3120];
% figure_range2 = [60:360];

z_distance   = 0.5;
ylim_range = [0 150];

temp_input_G1   = G1_Ca_binned(temp_range1,cell_id_G1);
temp_input_G2   = G2_Ca_binned(temp_range1,cell_id_G2);

% for training session
[~,num_cell_G1] = size(cell_id_G1);
[~,num_cell_G2] = size(cell_id_G2);

clear temp_input_fig_G1 % reset previous setting
for i = 1:num_cell_G1
    temp_input_fig_G1(:,i) = (temp_input_G1(:,i)) + i*z_distance;
end

clear temp_input_fig_G2 % reset previous setting
for i = 1:num_cell_G2
    temp_input_fig_G2(:,i) = (temp_input_G2(:,i)) + i*z_distance;
end

% for test session
% [~,num_cell] = size(cell_id);
% clear temp_input_fig % reset previous setting
% for i = 1:num_cell
%     temp_input_fig2(:,i) = (temp_input2(:,i)) + i*z_distance;
% end

%% ### for Fig.4c data export ###

% Use 1st (baseline), 3(US), 4(ITI-1), 5th(ITI-2) columns
for i = [1 3 4 5]
result_US_CA3_ctrl(:,i) = G1_results_resp0_removed_Z{i, 3}(1,:); % box plot for fig4d
result_US_CA3_ko(:,i)   = G2_results_resp0_removed_Z{i, 3}(1,:); % box plot for fig4d
end

%%


