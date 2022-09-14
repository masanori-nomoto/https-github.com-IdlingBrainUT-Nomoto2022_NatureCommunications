%% Cosine theta cal
% input time x neuron data
function [angle_raw, angle_mean] = fxn_nmt_cosine(reference_data, target_data)
% function [angle_raw, angle_mean, cos_theta_val, dot_val, bunbo] = fxn_nmt_cosine(reference_data, target_data) % for debug
%% for debug
% X = rand(10,7);
% Y = rand(5,7);

% X = [1 3];
% Y = [2 1]; % calculation results is 45deg.
% reference_data = X;
% target_data = Y;
%%
ref_mean = mean(reference_data,1);

target_data = mean(target_data,1); % target_mean on

dot_val = sum(ref_mean .* target_data,2);
%%
% cos_theta_val = dot_val/ (norm(ref_mean) .* norm(target_data,2));
bunbo = ( sqrt(sum((ref_mean.^2), 2)) .* sqrt(sum((target_data.^2),2)) );
cos_theta_val = dot_val./ bunbo;

cos_theta = acos(cos_theta_val); % radian
angle_raw = cos_theta *180/pi;
angle_mean = round(mean(angle_raw));
%%
% figure; histogram(angle_raw);
% hold on; histogram(angle_mean, 'FaceColor','r');
%% norm check
% norm_check = norm([1 -2 3])
% 
% test = (1^2) + (-2)^2 + (3^2);
% test2 = sqrt(test)
%%
end
%%