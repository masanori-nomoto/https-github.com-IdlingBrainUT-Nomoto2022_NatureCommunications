function ca_filt_data = fxn_filt_sigma_cutoff(ca_raw_data, sigma_cut_if_ON);

%% filtering code

    highpass = 0.01; % 30 sec highpass filter
    Ca_fs    = 20; % 20hz
    
    ca_temp_data = ca_raw_data   ;      %　vlockedデータの代入 ca 縦：時間　ｘ　横：細胞
   
        [b,a] = butter(1, highpass /(Ca_fs/2),'high'); %
        %[b,a] = cheby2(1, 10, highpass /(Ca_fs/2), 'high');% highpass =0.01 good  
        dataOut = filtfilt(b,a,ca_temp_data); %ゼロ位相filtfiltで ハイパス　フィルターを通す
        
%% zscoreも試す
% reply = input('Would you like to see an echo? (y/n): ','s');
% if strcmp(reply,'y')
%   disp(reply)
% end
% tf = strcmp(s1,s2) は、s1 と s2 を比較し、両者が同一の場合は 1 (true) を返し、そうでない場合は 0 (false) を返します。

% reply = input('Would you like to do Z-score normalization? (y/n): ','s');
reply = 'y'; % to automatically skip question
%%
if strcmp(reply,'y') % 
%   disp('You chose YES. Z-score normalization done!')
  
    dataIn1 = zscore(dataOut(1:12000,:)); 	%　縦にフィルターがかかる
    dataIn2 = zscore(dataOut(12001:48000,:)); 	%　縦にフィルターがかかる
    dataIn3 = zscore(dataOut(48001:55200,:)); 	%　縦にフィルターがかかる
%     dataIn4 = zscore(dataOut(55201:62400,:)); 	%　縦にフィルターがかかる
    
    dataIn = [dataIn1; dataIn2; dataIn3];   
    ca_filt_data = (dataIn);

    negative = find(ca_filt_data<sigma_cut_if_ON); % 3SD以下をカット
    ca_filt_data(negative) = zeros(size(negative));
%%  
elseif  strcmp(reply,'n') % 
  disp('You chose NO. Filtering and ZERO cutoff done!')
  
  ca_filt_data = (dataOut);

    negative = find(ca_filt_data<0); % 0
    ca_filt_data(negative) = zeros(size(negative));
else
   error('Unexpected situation')
end

%%
end