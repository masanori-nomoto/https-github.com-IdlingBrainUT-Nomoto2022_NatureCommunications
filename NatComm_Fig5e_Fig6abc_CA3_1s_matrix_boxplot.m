%% Responsiveness analysis
clc; clear; close all; disp('Previous data are cleared!')
%% load mat-file
load('NatComm_data_CA3_calcium_data_raw_ctrl11_ko9.mat');

%% generate concatenated data
% including trim data

G1_concatenated = [ca3_ctrl_m01, ca3_ctrl_m02, ca3_ctrl_m03, ca3_ctrl_m04, ca3_ctrl_m05, ...
                   ca3_ctrl_m06, ca3_ctrl_m07, ca3_ctrl_m08, ca3_ctrl_m09, ca3_ctrl_m10, ...
                   ca3_ctrl_m11]; 

G2_concatenated = [ca3_ko_m01, ca3_ko_m02, ca3_ko_m03, ca3_ko_m04, ca3_ko_m05, ...
                   ca3_ko_m06, ca3_ko_m07, ca3_ko_m08, ca3_ko_m09];

%% filtering ans sigma-cut off
% ca_filt_data = fxn_filt_sigma_ca(ca_raw_data, sigma_if_z_ON);

sigma_cut_if_ON = 3; % Z score cutoff 

G1_concatenated_temp  = fxn_ca3_filt_sigma_cutoff(G1_concatenated, sigma_cut_if_ON);
G2_concatenated_temp  = fxn_ca3_filt_sigma_cutoff(G2_concatenated, sigma_cut_if_ON);
%% Temporal binning
% [Ca_data_binned,  Ca_data_binned_mean] = funcHF_temporal_binning(Ca_data, binning_num)

binning_num = 20; % 20Hz to 1Hz data

[G1_Ca_binned,  G1_Ca_binned_mean] = fxn_temporal_binning(G1_concatenated_temp, binning_num);
[G2_Ca_binned,  G2_Ca_binned_mean] = fxn_temporal_binning(G2_concatenated_temp, binning_num);

reshape_range = (16:600);

% reshape_range = (151:600);
fxn_1s_bin_matrix(G1_Ca_binned, G2_Ca_binned, G1_Ca_binned_mean, G2_Ca_binned_mean, reshape_range)   
%% input time stamp information
upsample_rate = 1;

%% ##### Select all ITI-1 or narrow 9s ITI-1 (Latter 9s ITI-1 analysis)

% ### for normal analysis ###
% results = fxn_ca62400_12sessions(upsample_rate); % load session time information

% ### for iti1-latter 9s analysis ###
results = fxn_ca62400_12sessions_for_narrow_iti1(upsample_rate); % load session time information

%%
total_session_num  = 9;
responsiveness_threshold = 2;

G1_results = results;
G2_results = results;
disp('fxn_ca62400_12sessions done!')
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
        G1_results_responsive_mean  = cell2mat(G1_results_responsive(17:24,2:9)); % mean
        G1_results_responsive_error = cell2mat(G1_results_responsive(32:39,2:9)); % error
        
        G2_results_responsive_mean  = cell2mat(G2_results_responsive(17:24,2:9)); % mean
        G2_results_responsive_error = cell2mat(G2_results_responsive(32:39,2:9)); % error
        
        G1_G2_results_responsive_subtracttion = G1_results_responsive_mean-G2_results_responsive_mean;
%% mean, error, and ranksum-stat figure
figure;
subplot(232); imagesc(G1_results_responsive_mean)
subplot(233); imagesc(G2_results_responsive_mean)
subplot(234); imagesc(G1_G2_results_responsive_subtracttion)
subplot(235); imagesc(G1_G2_results_responsive_ranksum_p)
subplot(236); imagesc(G1_G2_results_responsive_ranksum_h)
%% Baseline-responsiveness zero-cell removed Reponsive 1s bin matrix (ALL)
% (ALL cells resutls, removed IDs showeing infinite value, due to baseline zero value)

reshape_range   = (16:600);
% reshape_range = (151:600);

G1_Ca_binned_responsive = G1_results{1,6}(:,G1_results{2,7}); G1_Ca_binned_responsive_mean = mean(G1_Ca_binned_responsive,2);
G2_Ca_binned_responsive = G2_results{1,6}(:,G2_results{2,7}); G2_Ca_binned_responsive_mean = mean(G2_Ca_binned_responsive,2);

fxn_1s_bin_matrix(G1_Ca_binned_responsive, G2_Ca_binned_responsive, ...
                      G1_Ca_binned_responsive_mean, G2_Ca_binned_responsive_mean, reshape_range)
%% Baseline-responsiveness zero-cell removed Reponsive 1s bin matrix (CS,US, and LTM-CS)
% (CS,US, and LTM-CS cells resutls, removed IDs showeing infinite value, due to baseline zero value)
% CS, US, LTM-CS cell pickup

reshape_range   = (16:600);
% reshape_range = (151:600);

G1_neruon_ID = unique([G1_results{2,5}, G1_results{3,5}, G1_results{11,5}]);
G2_neruon_ID = unique([G2_results{2,5}, G2_results{3,5}, G2_results{11,5}]);

G1_Ca_binned_responsive = G1_results{1,6}(:,G1_neruon_ID); G1_Ca_binned_responsive_mean = mean(G1_Ca_binned_responsive,2);
G2_Ca_binned_responsive = G2_results{1,6}(:,G2_neruon_ID); G2_Ca_binned_responsive_mean = mean(G2_Ca_binned_responsive,2);

fxn_1s_bin_matrix(G1_Ca_binned_responsive, G2_Ca_binned_responsive, ...
                      G1_Ca_binned_responsive_mean, G2_Ca_binned_responsive_mean, reshape_range)   
%% Box plot
% fxn_boxplot_double(G1_results_ratio, G2_results_ratio, input_session_num, section_num);
section_num_end = 9; % con-ko -> 14

for i = 2:section_num_end
 fxn_boxplot_double(G1_results_responsive, G2_results_responsive, total_session_num, i);
end
%%
close all
% ##### Z-score data using 

%% Baseline-responsiveness zero-cell removed Z-score 1s bin matrix
% (ALL cells resutls, removed IDs showeing infinite value, due to baseline zero value)

reshape_range   = (16:600);
% reshape_range = (151:600);

G1_Ca_binned_resp0_removed_Z = G1_Ca_binned(:,G1_results{2,7}); G1_Ca_binned_resp0_removed_Z_mean = mean(G1_Ca_binned_resp0_removed_Z,2);
G2_Ca_binned_resp0_removed_Z = G2_Ca_binned(:,G2_results{2,7}); G2_Ca_binned_resp0_removed_Z_mean = mean(G2_Ca_binned_resp0_removed_Z,2);

fxn_1s_bin_matrix(G1_Ca_binned_resp0_removed_Z, G2_Ca_binned_resp0_removed_Z, ...
                      G1_Ca_binned_resp0_removed_Z_mean, G2_Ca_binned_resp0_removed_Z_mean, reshape_range)
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
            G1_G2_results_resp0_removed_Z_ranksum_stats_zval(ii,i) = G1_G2_results_resp0_removed_Z_ranksum_stats(ii, i).zval    ;
    end
end
        G1_results_resp0_removed_Z_mean  = cell2mat(G1_results_resp0_removed_Z(16:24,2:9)); % mean
        G1_results_resp0_removed_Z_error = cell2mat(G1_results_resp0_removed_Z(31:39,2:9)); % error
        
        G2_results_resp0_removed_Z_mean  = cell2mat(G2_results_resp0_removed_Z(16:24,2:9)); % mean
        G2_results_resp0_removed_Z_error = cell2mat(G2_results_resp0_removed_Z(31:39,2:9)); % error
        
        G1_G2_results_resp0_removed_Z_subtracttion = G1_results_resp0_removed_Z_mean-G2_results_resp0_removed_Z_mean;
        %% Signrank calculaiton with removed Z-score matrix 
for i = 2:total_session_num
    for ii = 1:total_session_num
        for iii = 1:1:total_session_num
% [p,h,stats] = signrank(___)
[G1_results_resp0_removed_Z_signrank_p{iii,ii,i}, G1_results_resp0_removed_Z_signrank_h{iii,ii,i} ...
G1_results_resp0_removed_Z_signrank_stats{iii,ii,i}] = signrank((G1_results_resp0_removed_Z{ii,i}),(G1_results_resp0_removed_Z{iii,i})); 

[G2_results_resp0_removed_Z_signrank_p{iii,ii,i}, G2_results_resp0_removed_Z_signrank_h{iii,ii,i} ...
G2_results_resp0_removed_Z_signrank_stats{iii,ii,i}] = signrank((G2_results_resp0_removed_Z{ii,i}),(G2_results_resp0_removed_Z{iii,i})); 
        
        end
    end
end
%% export signrank stats table

G1_results_resp0_removed_Z_signrank_CS_p     = G1_results_resp0_removed_Z_signrank_p(:,:,2);
G1_results_resp0_removed_Z_signrank_CS_stats = G1_results_resp0_removed_Z_signrank_stats(:,:,2);

G1_results_resp0_removed_Z_signrank_US_p     = G1_results_resp0_removed_Z_signrank_p(:,:,3);
G1_results_resp0_removed_Z_signrank_US_stats = G1_results_resp0_removed_Z_signrank_stats(:,:,3);

G1_results_resp0_removed_Z_signrank_LTMCS_p     = G1_results_resp0_removed_Z_signrank_p(:,:,8);
G1_results_resp0_removed_Z_signrank_LTMCS_stats = G1_results_resp0_removed_Z_signrank_stats(:,:,8);

for i = 2:total_session_num
    for ii = 1:total_session_num
        try    
    G1_results_resp0_removed_Z_signrank_CS_zval(ii,i) = G1_results_resp0_removed_Z_signrank_CS_stats{ii, i}.zval ;
    G1_results_resp0_removed_Z_signrank_US_zval(ii,i) = G1_results_resp0_removed_Z_signrank_US_stats{ii, i}.zval ;
    G1_results_resp0_removed_Z_signrank_LTMCS_zval(ii,i) = G1_results_resp0_removed_Z_signrank_LTMCS_stats{ii, i}.zval ;
        catch
        end
    end
end

G2_results_resp0_removed_Z_signrank_CS_p     = G2_results_resp0_removed_Z_signrank_p(:,:,2);
G2_results_resp0_removed_Z_signrank_CS_stats = G2_results_resp0_removed_Z_signrank_stats(:,:,2);

G2_results_resp0_removed_Z_signrank_US_p     = G2_results_resp0_removed_Z_signrank_p(:,:,3);
G2_results_resp0_removed_Z_signrank_US_stats = G2_results_resp0_removed_Z_signrank_stats(:,:,3);

G2_results_resp0_removed_Z_signrank_LTMCS_p     = G2_results_resp0_removed_Z_signrank_p(:,:,8);
G2_results_resp0_removed_Z_signrank_LTMCS_stats = G2_results_resp0_removed_Z_signrank_stats(:,:,8);

for i = 2:total_session_num
    for ii = 1:total_session_num
        try    
    G2_results_resp0_removed_Z_signrank_CS_zval(ii,i) = G2_results_resp0_removed_Z_signrank_CS_stats{ii, i}.zval ;
    G2_results_resp0_removed_Z_signrank_US_zval(ii,i) = G2_results_resp0_removed_Z_signrank_US_stats{ii, i}.zval ;
    G2_results_resp0_removed_Z_signrank_LTMCS_zval(ii,i) = G2_results_resp0_removed_Z_signrank_LTMCS_stats{ii, i}.zval ;
        catch
        end
    end
end
%% Baseline-responsiveness zero-cell removed Z-score 1s bin matrix  
% (CS,US, and LTM-CS cells resutls, removed IDs showeing infinite value, due to baseline zero value)
% CS, US, LTM-CS cell pickup

reshape_range   = (16:600);
% reshape_range = (151:600);

% G1_neruon_ID = unique([G1_results{2,5}]); % CS
% G2_neruon_ID = unique([G2_results{2,5}]); % CS

% G1_neruon_ID = unique([G1_results{3,5}]); % US
% G2_neruon_ID = unique([G2_results{3,5}]); % US

G1_neruon_ID = unique([G1_results{8,5}]); % LTM-CS
G2_neruon_ID = unique([G2_results{8,5}]); % LTM-CS

% G1_neruon_ID = unique([G1_results{2,5}, G1_results{3,5}, G1_results{8,5}]); % CS,US, and LTM-CS
% G2_neruon_ID = unique([G2_results{2,5}, G2_results{3,5}, G2_results{8,5}]); % CS,US, and LTM-CS

G1_Ca_binned_resp0_removed_Z = G1_Ca_binned(:,G1_neruon_ID); G1_Ca_binned_resp0_removed_Z_mean = mean(G1_Ca_binned_resp0_removed_Z,2);
G2_Ca_binned_resp0_removed_Z = G2_Ca_binned(:,G2_neruon_ID); G2_Ca_binned_resp0_removed_Z_mean = mean(G2_Ca_binned_resp0_removed_Z,2);

fxn_1s_bin_matrix(G1_Ca_binned_resp0_removed_Z, G2_Ca_binned_resp0_removed_Z, ...
                      G1_Ca_binned_resp0_removed_Z_mean, G2_Ca_binned_resp0_removed_Z_mean, reshape_range) 
%% Box plot
% fxn_boxplot_double(G1_results_ratio, G2_results_ratio, input_session_num, section_num);
section_num_end = 9; % 9=STM removed,  12=Rest summarised, 14=entire

for i = 1:section_num_end
 fxn_boxplot_double(G1_results_resp0_removed_Z, G2_results_resp0_removed_Z, section_num_end, i);
% ylim([-0.05 0.6]); % for CS and test-cs
ylim([-0.05 1.2]); % for CS and test-cs
box off; 
end

%% ### LTM-test activity export for figure 5e

% for ctrl
G1_LTM_test_CS_acc_values(:,1) = (G1_results_resp0_removed_Z{7, 2});
G1_LTM_test_CS_acc_values(:,2) = (G1_results_resp0_removed_Z{8, 2});
G1_LTM_test_CS_acc_values(:,3) = (G1_results_resp0_removed_Z{9, 2});

G1_LTM_test_US_acc_values(:,1) = (G1_results_resp0_removed_Z{7, 3});
G1_LTM_test_US_acc_values(:,2) = (G1_results_resp0_removed_Z{8, 3});
G1_LTM_test_US_acc_values(:,3) = (G1_results_resp0_removed_Z{9, 3});

G1_LTM_test_LTMCS_acc_values(:,1) = (G1_results_resp0_removed_Z{7, 8});
G1_LTM_test_LTMCS_acc_values(:,2) = (G1_results_resp0_removed_Z{8, 8});
G1_LTM_test_LTMCS_acc_values(:,3) = (G1_results_resp0_removed_Z{9, 8});

% for ko
G2_LTM_test_CS_acc_values(:,1) = (G2_results_resp0_removed_Z{7, 2});
G2_LTM_test_CS_acc_values(:,2) = (G2_results_resp0_removed_Z{8, 2});
G2_LTM_test_CS_acc_values(:,3) = (G2_results_resp0_removed_Z{9, 2});

G2_LTM_test_US_acc_values(:,1) = (G2_results_resp0_removed_Z{7, 3});
G2_LTM_test_US_acc_values(:,2) = (G2_results_resp0_removed_Z{8, 3});
G2_LTM_test_US_acc_values(:,3) = (G2_results_resp0_removed_Z{9, 3});

G2_LTM_test_LTMCS_acc_values(:,1) = (G2_results_resp0_removed_Z{7, 8});
G2_LTM_test_LTMCS_acc_values(:,2) = (G2_results_resp0_removed_Z{8, 8});
G2_LTM_test_LTMCS_acc_values(:,3) = (G2_results_resp0_removed_Z{9, 8});

%% ### LTM-test activity export for figure 6abc.

% for G1_results, for ctrl
for i = 1:6
G1_boxplot_6sessions_CS_cells_values(:,i) = (G1_results_resp0_removed_Z{i, 2});
G1_boxplot_6sessions_US_cells_values(:,i) = (G1_results_resp0_removed_Z{i, 3});
G1_boxplot_6sessions_LTMCS_cells_values(:,i) = (G1_results_resp0_removed_Z{i, 8});
end

% for G2_results, for ko
for i = 1:6
G2_boxplot_6sessions_CS_cells_values(:,i) = (G2_results_resp0_removed_Z{i, 2});
G2_boxplot_6sessions_US_cells_values(:,i) = (G2_results_resp0_removed_Z{i, 3});
G2_boxplot_6sessions_LTMCS_cells_values(:,i) = (G2_results_resp0_removed_Z{i, 8});
end
%%