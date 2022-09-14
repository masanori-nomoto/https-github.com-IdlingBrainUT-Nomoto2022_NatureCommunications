%% Maharanobis Population vector distance calculation code
% 1bin_ca_data: time x neuron matrix
function [res_MPPCA, res_thrcov_PCA] = fxn_Marchenko_thrcover_PCA(ca_data, thrcov_PC_percent)
%% for debug,
% PC_threshold = 90; % 90 means PCs compose 90 percentile 
X = ca_data;
PC_threshold = thrcov_PC_percent;
%% MP-PCA dimension reduction
% n × B matrix, (1+sqrt(n/B)).^2
% PCA calculation, [coeff,score,latent] = pca(X);

% Marchenko–Pastur cal
[size_X(1,1), size_X(1,2)] = size(X); % dim1 -> time, dim2 -> neuron
% [size_Y(1,1), size_Y(1,2)] = size(Y); % dim1 -> time, dim2 -> neuron

% (1+sqrt(n/t)).^2; % n -> neuron, t -> time;
[PCA_X_coeff, PCA_X_score, PCA_X_latent] = pca(X);
n1 = size_X(1,2); t1 = size_X(1,1);
MP_val_X = (1+sqrt(n1/t1)).^2; % n -> neuron, t -> time;
PCA_X_sig = find(PCA_X_latent > MP_val_X);
PCA_X_MPsig_num = numel(PCA_X_sig);

PCA_X_latent_cumsum(:,1) = PCA_X_latent;
PCA_X_latent_cumsum(:,2) = cumsum(PCA_X_latent);
PCA_X_latent_cumsum(:,3) = PCA_X_latent_cumsum(:,2)/sum(PCA_X_latent)*100;
PCA_X_latent_cumsum(:,4) = PCA_X_latent_cumsum(:,3) < PC_threshold;

% PCA_X_PC_thr_num = sum(PCA_X_latent_cumsum(:,4),1); % ### original ###
PCA_X_PC_thr_num = 3; % #### forced top3 PCs ON ###

MPPCA_X_score     = PCA_X_score(:,1:PCA_X_MPsig_num);
thr_cover_X_score = PCA_X_score(:,1:PCA_X_PC_thr_num);

% [PCA_Y_coeff, PCA_Y_score, PCA_Y_latent] = pca(Y);
% n2 = size_Y(1,2); t2 = size_Y(1,1);
% MP_val_Y = (1+sqrt(n2/t2)).^2; % n -> neuron, t -> time;
% PCA_Y_sig = find(PCA_Y_latent > MP_val_Y);
% PCA_Y_MPsig_num = numel(PCA_Y_sig);
% PCA_Y_latent_cumsum(:,1) = PCA_Y_latent;
% PCA_Y_latent_cumsum(:,2) = cumsum(PCA_Y_latent);
% PCA_Y_latent_cumsum(:,3) = PCA_Y_latent_cumsum(:,2)/sum(PCA_Y_latent)*100;
% PCA_Y_latent_cumsum(:,4) = PCA_Y_latent_cumsum(:,3) < PC_threshold;
% PCA_Y_PC_thr_num = sum(PCA_Y_latent_cumsum(:,4),1);

% figure('Position',[600,50,400,200]); 
% % subplot(211); 
% plot(PCA_X_latent,'k'); hold on
% area(PCA_X_MPsig_num, max(PCA_X_latent),'EdgeColor', 'red','FaceColor', 'red'); 
% title('Marchenko–Pastur Max Lambda for reference matrix'); ylabel('Eigenvalue'); xlabel('Latent')
% area(PCA_X_PC_thr_num, max(PCA_X_latent),'EdgeColor','blue','FaceColor', 'blue'); 
% title('Dimension reduction thresholding'); ylabel('Eigenvalue'); xlabel('Latent')
% legend PCA-latent-eigenval  Marchenko–Pastur-thr PC-coverage-thr
% 
% ax = gca;
% set(gca, 'FontSize', 10, 'FontName','Arial'); colormap('parula'); grid on; ax.TickDir = 'both';

% subplot(212); plot(PCA_Y_latent,'k'); hold on
% area(PCA_Y_MPsig_num, max(PCA_Y_latent),'EdgeColor', 'red','FaceColor', 'red'); 
% title('Marchenko–Pastur Max Lambda for reference matrix'); ylabel('Eigenvalue'); xlabel('Latent')
% area(PCA_Y_PC_thr_num, max(PCA_Y_latent),'EdgeColor','blue','FaceColor', 'blue'); 
% title('Marchenko–Pastur Max Lambda for reference matrix'); ylabel('Eigenvalue'); xlabel('Latent')
% legend PCA-latent  Marchenko–Pastur-sig PC-coverage

% hold off
%%

%% 
res_MPPCA.MPPCA_sig_num          = PCA_X_MPsig_num;
res_MPPCA.MPPCA_sig_score           = MPPCA_X_score;
res_MPPCA.thrcov_PCA_latent_cumsum     = PCA_X_latent_cumsum;
res_MPPCA.MP_Lambda_max_val     = MP_val_X;

res_thrcov_PCA.thrcov_PCA_thr_num        = PCA_X_PC_thr_num;
res_thrcov_PCA.thrcov_PCA_score          = thr_cover_X_score;
res_thrcov_PCA.thrcov_PCA_latent_cumsum  = PCA_X_latent_cumsum;
res_thrcov_PCA.thrcov_PCA_threshold      =  PC_threshold;
%%
end
%%
