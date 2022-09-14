function [results] = fxnHF_ca57600_8sessions(onset_bin)

%% session information in second scale
%        res_temp     = ToIntervals(time_stamp, A1(:,i_A1));
%         A1_table = cat(1,A1_table,res_temp);
% onset_bin1 = 360;

% t.s1 = [1 :onset_bin+150]; % baseline
t.s1 = [onset_bin+1 :onset_bin+150]; % baseline
t.s2 = []; % CS 
t.s3 = []; % US
t.s4 = []; % ITI-1
t.s5 = []; % ITI-2
% t.s6 = [601:2400]; % Homecage-1
% t.s7 = [2401:2520]; % STM baseline
% t.s8 = [2521:2550 2581:2610 2641:2670 2701:2730]; % STM CS
% t.s9 = [2551:2580 2611:2640 2671:2700 2731:2760]; % STM ITI
t.s6 = [2521:2640]; % LTM baseline
t.s7 = [2641:2670 2701:2730 2761:2790 2821:2850]; % LTM CS
t.s8 = [2671:2700 2731:2760 2791:2820 2851:2880]; % LTM ITI

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

for i = 1:10
t_temp_s2 = [onset_bin+151+45*(i-1): onset_bin+151+45*(i-1)+14]; t.s2 = cat(2,t.s2,t_temp_s2); % for CS only
% t_temp_s2 = [onset_bin+151+45*(i-1): onset_bin+151+45*(i-1)+11]; t.s2 = cat(2,t.s2,t_temp_s2); % for CS-US
t_temp_s3 = [onset_bin+163+45*(i-1): onset_bin+163+45*(i-1)+2];  t.s3 = cat(2,t.s3,t_temp_s3); % US
t_temp_s4 = [onset_bin+166+45*(i-1): onset_bin+166+45*(i-1)+14]; t.s4 = cat(2,t.s4,t_temp_s4); % ITI-E
t_temp_s5 = [onset_bin+181+45*(i-1): onset_bin+181+45*(i-1)+14]; t.s5 = cat(2,t.s5,t_temp_s5); % ITI-L
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
% results{9,1} = t.s9;
% results{10,1} = t.s10;
% results{11,1} = t.s11;
% results{12,1} = t.s12;
% results{13,1} = t.s13;
% results{14,1} = t.s14;

%%
end