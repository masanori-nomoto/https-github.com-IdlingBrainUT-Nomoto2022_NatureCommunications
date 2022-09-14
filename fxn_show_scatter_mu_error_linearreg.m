%% Show scatter mean figure
function fxn_show_scatter_mu_error_linearreg(fig_xy_dim ,data1_all_raw, data2_all_raw, xlabel_title_input, ylabel_title_input) 

%% for debug
% demo_data1 = rand(10,2).*2;
% demo_data2 = rand(20,2).*10;
% fig_xy_dim = [140,130] % def [140 100]
%% 
demo_data1 = data1_all_raw;
demo_data2 = data2_all_raw;
%%
[data1_mu, data1_error] = fxn_mu_error(demo_data1);
[data2_mu, data2_error] = fxn_mu_error(demo_data2);
%%
mu = [data1_mu; data2_mu];
error = [data1_error; data2_error];
%%
data1_ave_color = '#0000ff'; % blue
data2_ave_color = '#ff0000'; % red

data1_all_color = '#4682d4'; % gray blue
data1_all_color_edge = '#192F4D'; % dark blue

data2_all_color = '#db7093'; % gray red
data2_all_color_edge = '#800D33'; % dark red

figure('Position',[200,150, fig_xy_dim]); %[left bottom width height]
%%

er1 = errorbar(mu(1,1),mu(1,2),error(1,2),error(1,2),error(1,1),error(1,1), ...
    'LineWidth', 1, 'MarkerSize',2, 'Marker','o','Color', data1_ave_color, ...
    'MarkerFaceColor',	data1_ave_color, 'MarkerEdgeColor',data1_ave_color);

% Set transparency level (0:1)
alpha = 0.65;   
% Set transparency (undocumented)
set([er1.Bar, er1.Line], 'ColorType', 'truecoloralpha', 'ColorData', [er1.Line.ColorData(1:3); 255*alpha])

hold on

h = plot(demo_data1(:,1), demo_data1(:,2), '.', ...
    'LineWidth', 0.5, 'Color', data1_all_color, 'MarkerSize',8, ...
    'MarkerFaceColor',	data1_all_color_edge, 'MarkerEdgeColor',data1_all_color); 

hold on;

x1 = demo_data1(:,1); y1 = demo_data1(:,2);
axis([min(x1) max(x1) min(y1) max(y1)]); % desired limit for line
lsline; 
% xlim auto% This is equivalent to refline(ax1)
%%

hold on

er2 = errorbar(mu(2,1),mu(2,2),error(2,2),error(2,2),error(2,1),error(2,1), ...
    'LineWidth', 1, 'MarkerSize',2, 'Marker','o','Color', data2_ave_color, ...
    'MarkerFaceColor',	data2_ave_color, 'MarkerEdgeColor',data2_ave_color); 

% Set transparency level (0:1)
alpha = 0.65;   
% Set transparency (undocumented)
set([er2.Bar, er2.Line], 'ColorType', 'truecoloralpha', 'ColorData', [er2.Line.ColorData(1:3); 255*alpha])

hold on;

plot(demo_data2(:,1), demo_data2(:,2), '.', ...
    'LineWidth', 0.5, 'Color', data2_all_color, 'MarkerSize',8, ...
    'MarkerFaceColor',	data2_all_color_edge, 'MarkerEdgeColor',data2_all_color); 

x2 = demo_data2(:,1); y2 = demo_data2(:,2);
axis([min(x2) max(x2) min(y2) max(y2)]); % desired limit for line
set(h, 'HandleVisibility', 'off'); % hide previous plot
lsline % plot line
set(h, 'HandleVisibility', 'on'); % restore visibility
xlim auto % restore axis limit

%% ticks control
% for debug
% xtickname1 = 'CS'
% xtickname2 = 'ITI'
% ylabel_title_input =  {'Relative similality of';'PVD to ITI'};

xlabel(xlabel_title_input); 
ylabel(ylabel_title_input);

xlim([0 max([x1; x2])]);
ylim([0 max([y1; y2])]);
% yticks([0:10:60]);

box off
set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 7, 'FontName','Arial'); %grid on;

box off
%%

end
%%
