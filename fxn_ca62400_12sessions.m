function [results] = fxn_ca62400_12sessions(upsample_rate)
%% for debug,
% upsample 1s time stamp num -> 4 x 50ms = 200ms binning
% upsample 1s time stamp num -> 1 = 1 s binning
% upsample_rate = 1; % upsample 1s time stamp num -> 50ms x 4 = 200ms binning

%% session information in second scale
%        res_temp     = ToIntervals(time_stamp, A1(:,i_A1));
%         A1_table = cat(1,A1_table,res_temp);


t.s1 = [(0*upsample_rate)+1:150*upsample_rate]; % baseline
t.s2 = []; % CS 
t.s3 = []; % US
t.s4 = []; % ITI-1
t.s5 = []; % ITI-2
t.s6 =  [600*upsample_rate+1:2400*upsample_rate]; % Homecage-1
t.s7 =  [2400*upsample_rate+1:2520*upsample_rate]; % STM baseline
t.s8 =  [2520*upsample_rate+1:2550*upsample_rate 2580*upsample_rate+1:2610*upsample_rate 2640*upsample_rate+1:2670*upsample_rate 2700*upsample_rate+1:2730*upsample_rate]; % STM CS
t.s9 =  [2550*upsample_rate+1:2580*upsample_rate 2610*upsample_rate+1:2640*upsample_rate 2670*upsample_rate+1:2700*upsample_rate 2730*upsample_rate+1:2760*upsample_rate]; % STM ITI
t.s10 = [2760*upsample_rate+1:2880*upsample_rate]; % LTM baseline
t.s11 = [2880*upsample_rate+1:2910*upsample_rate 2940*upsample_rate+1:2970*upsample_rate 3000*upsample_rate+1:3030*upsample_rate 3060*upsample_rate+1:3090*upsample_rate]; % LTM CS
t.s12 = [2910*upsample_rate+1:2940*upsample_rate 2970*upsample_rate+1:3000*upsample_rate 3030*upsample_rate+1:3060*upsample_rate 3090*upsample_rate+1:3120*upsample_rate]; % LTM ITI

% t.s1 = [1:150]; % baseline
% t.s2 = []; % CS 
% t.s3 = []; % US
% t.s4 = []; % ITI-1
% t.s5 = []; % ITI-2
% t.s6 = [601:2400]; % Homecage-1
% % t.s7 = [1201:1800]; % Homecage-2
% % t.s8 = [1801:2400]; % Homecage-3
% t.s9 = [2401:2520]; % STM baseline
% t.s10 = [2521:2550 2581:2610 2641:2670 2701:2730]; % STM CS
% t.s11 = [2551:2580 2611:2640 2671:2700 2731:2760]; % STM ITI
% t.s12 = [2761:2880]; % LTM baseline
% t.s13 = [2881:2910 2941:2970 3001:3030 3061:3090]; % LTM CS
% t.s14 = [2911:2940 2971:3000 3031:3060 3091:3120]; % LTM ITI

for i = 1:10 % default
% t_temp_s2 = [151+45*(i-1): 151+45*(i-1)+11]; t.s2 = cat(2,t.s2,t_temp_s2); % CS
% t_temp_s3 = [163+45*(i-1): 163+45*(i-1)+2];  t.s3 = cat(2,t.s3,t_temp_s3); % US
% t_temp_s4 = [166+45*(i-1): 166+45*(i-1)+14]; t.s4 = cat(2,t.s4,t_temp_s4); % ITI-E
% t_temp_s5 = [181+45*(i-1): 181+45*(i-1)+14]; t.s5 = cat(2,t.s5,t_temp_s5); % ITI-L
t_temp_s2 = [150*upsample_rate+1 + 45*upsample_rate*(i-1): 150*upsample_rate + 45*upsample_rate*(i-1) + 12*upsample_rate]; 
t.s2 = cat(2,t.s2,t_temp_s2); % CS
t_temp_s3 = [162*upsample_rate+1 + 45*upsample_rate*(i-1): 162*upsample_rate + 45*upsample_rate*(i-1) + 3*upsample_rate];
t.s3 = cat(2,t.s3,t_temp_s3); % US
t_temp_s4 = [165*upsample_rate+1 + 45*upsample_rate*(i-1): 165*upsample_rate + 45*upsample_rate*(i-1) + 15*upsample_rate];
t.s4 = cat(2,t.s4,t_temp_s4); % ITI-1
t_temp_s5 = [180*upsample_rate+1 + 45*upsample_rate*(i-1): 180*upsample_rate + 45*upsample_rate*(i-1) + 15*upsample_rate]; 
t.s5 = cat(2,t.s5,t_temp_s5); % ITI-2
end

%% frame sorting
results={}; % initialize

results{1,1} = t.s1;
results{2,1} = t.s2;
results{3,1} = t.s3;
results{4,1} = t.s4;
results{5,1} = t.s5;
results{6,1} = t.s6;
results{7,1} = t.s7;
results{8,1} = t.s8;
results{9,1} = t.s9;
results{10,1} = t.s10;
results{11,1} = t.s11;
results{12,1} = t.s12;
% results{13,1} = t.s13;
% results{14,1} = t.s14;

%%
end