%% Calculate Maharanobis distance
% this code calculate Maharanobis distance by using covariance matrix of reference data
% input time x neuron matrix
function [Mahalanobis_d_raw, Mahalanobis_mean] = fxn_nmt_Mahalanobis(reference_data, target_data)
%% for debug,
% input_data = score_MPPCA;
% reference_data = input_data(100:300,:); 
% target_data    = input_data([100:300 1000:3000],:);
%%
X = reference_data;
Y = target_data;
%% debug
% figure;
% subplot(211); plot(X);
% subplot(212); plot(Y);
%% Calculate Maharanobis distance
%  M^2 = (X(I,:)-MU) * SIGMA^(-1) * (X(I,:)-MU)',
% X = time x neuron, reference
% Y = time x neuron, target

muX = mean(X,1); % n dimention

% covariance_matrix type
covariance_X = cov(X);
% covariance_X = cov([X;Y]);

inv_covX = inv(covariance_X);

% % b*inv(A) を b/A に置き換えます。
M2 = (Y-muX) * inv_covX * (Y-muX)'; % M2 not squared
M2_diag = diag(M2);
Mahalanobis_d_raw = sqrt(M2_diag);
Mahalanobis_mean = mean(Mahalanobis_d_raw,1);
%%
% openExample('stats/CompareMahalanobisAndSquaredEuclideanDistancesExample')
% Mahal_demo = mahal(Y,X);
%%
%%
% figure;plot(Mahalanobis_d_raw*100); hold on
% plot(Mahal_demo);

end