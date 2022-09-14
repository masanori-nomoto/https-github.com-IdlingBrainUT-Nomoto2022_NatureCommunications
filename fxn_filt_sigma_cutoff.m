function ca_filt_data = fxn_filt_sigma_cutoff(ca_raw_data, sigma_cut_if_ON);

%% filtering code

    highpass = 0.01; % highpass filter
    Ca_fs    = 20; % 20hz
    
    ca_temp_data = ca_raw_data   ;      %　time x cells
   
        [b,a] = butter(1, highpass /(Ca_fs/2),'high'); %
        %[b,a] = cheby2(1, 10, highpass /(Ca_fs/2), 'high');% highpass =0.01 good  
        dataOut = filtfilt(b,a,ca_temp_data); %zero filtfilt
        

reply = 'y'; % skipping
%%
if strcmp(reply,'y') % 
  disp('You chose YES. Z-score normalization done!')
  
    dataIn1 = zscore(dataOut(1:12000,:)); 	%　
    dataIn2 = zscore(dataOut(12001:48000,:)); 	%　
    dataIn3 = zscore(dataOut(48001:55200,:)); 	%　
    dataIn4 = zscore(dataOut(55201:62400,:)); 	%　
    
    dataIn = [dataIn1; dataIn2; dataIn3; dataIn4];   
    ca_filt_data = (dataIn);

    negative = find(ca_filt_data<sigma_cut_if_ON); % 3SD cutoff
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