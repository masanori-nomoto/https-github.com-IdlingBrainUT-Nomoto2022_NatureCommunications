%% figure line plot
function [ttest_results] = fxn_line_plot_for_imported_data(data1_all_raw, data2_all_raw, xtickname1, xtickname2, ylabel_title_input) 
% input samples x time point data

%% for debug
% xtickname1 = 'CS'
% xtickname2 = 'ITI'
% ylabel_title_input =  {'Relative similality of';'PVD to ITI'};
% data1_all = mahal_ca1_cs_iti_all_wt'    ; % samples x time point
% data2_all = mahal_ca1_cs_iti_all_ko'    ; % samples x time point
%% data import
% input samples x time point
data1_all = data1_all_raw'    ; % samples x time point
data2_all = data2_all_raw'    ; % samples x time point
%% data process
data1_all_ave = mean(data1_all,2);
data2_all_ave = mean(data2_all,2);

w = 0;
data1_all_error = std(data1_all,w,2)./sqrt(size(data1_all,2));
data2_all_error = std(data2_all,w,2)./sqrt(size(data2_all,2));

%% figure data

data1_ave_color = '#0000ff'; % blue
data2_ave_color = '#ff0000'; % red

data1_all_color = '#b0c4de'; % gray blue
data2_all_color = '#ffc0cb'; % gray red

figure('Position',[200,150, 70,100]); %[left bottom width height]

% first individual
plot(data1_all(:,:), 'LineWidth', 0.4, ...
    'Marker','.','Color',data1_all_color,'MarkerFaceColor',	data1_all_color, 'MarkerEdgeColor',data1_all_color)
hold on
plot(data2_all(:,:),'LineWidth', 0.4, ... 
    'Marker','.','Color',data2_all_color,'MarkerFaceColor',	data2_all_color, 'MarkerEdgeColor',data2_all_color)
hold on

% second average
errorbar(data1_all_ave(:,:), data1_all_error,'LineWidth', 1.2, 'MarkerSize',3, ...
    'Marker','o','Color',data1_ave_color,'MarkerFaceColor',	data1_ave_color, 'MarkerEdgeColor',data1_ave_color);
hold on
errorbar(data2_all_ave(:,:), data2_all_error, 'LineWidth', 1.2, 'MarkerSize',3, ... 
    'Marker','o','Color',data2_ave_color,'MarkerFaceColor',	data2_ave_color, 'MarkerEdgeColor',data2_ave_color);
hold off
% ticks control
xlim([0.7 2.3]) ;xticks([1:3]); xticklabels({xtickname1,xtickname2});

% yticks([0:0.2:1.5]);
ylabel(ylabel_title_input); 
box off
set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 7, 'FontName','Arial'); %grid on;

%% ttest

x = data1_all_raw(:,2);
y = data2_all_raw(:,2);

[h,p,ci,stats] = ttest2(x,y);

 ttest_results.h     = h;
 ttest_results.p     = p;
 ttest_results.ci    = ci;
 ttest_results.stats = stats; 
%%

end

