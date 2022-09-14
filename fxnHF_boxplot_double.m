function fxnHF_boxplot_double(G1_results_ratio, G2_results_ratio, input_session_num, section_num,box_name_cell);
%%
if input_session_num == 4
% section_num = 2;
Box_X = [];
Box_Y = [];


for i = 1:input_session_num
        cat_x_temp_G1 = (cell2mat(G1_results_ratio(i,section_num)))';
        cat_x_temp_G2 = (cell2mat(G2_results_ratio(i,section_num)))';
    Box_X = [Box_X; cat_x_temp_G1];
    Box_X = [Box_X; cat_x_temp_G2];  
    
             if i == 1; group1 = ('Baseline');  group2 = ('Baseline-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 2; group1 = ('CS');  group2 = ('CS-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 3; group1 = ('US');  group2 = ('US-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];   
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 4; group1 = ('ITI-1');  group2 = ('ITI-1-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 5; group1 = ('ITI-2');  group2 = ('ITI-2-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 6; group1 = ('Rest-1');  group2 = ('Rest-1-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 7; group1 = ('HC-2');  group2 = ('HC-2-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 8; group1 = ('HC-3');  group2 = ('HC-3-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 7; group1 = ('STM-Baseline');  group2 = ('STM-Baseline-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 8; group1 = ('STM-CS');  group2 = ('STM-CS-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 9; group1 = ('STM-ITI');  group2 = ('STM-ITI-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 6; group1 = ('LTM-Baseline');  group2 = ('LTM-baseline-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 7; group1 = ('LTM-CS'); group2 = ('LTM-CS-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 8; group1 = ('LTM-ITI'); group2 = ('LTM-ITI-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];    
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         end  
    end

%% Boxplot new
% figure('color',[1,1,1]);
figure('Position',[244,181,160,130]); %[left bottom width height]
boxplot(Box_X, Box_Y,'factorgap',2,'color','kr','OutlierSize',2, 'Widths' ,1.0,'Symbol','','BoxStyle','outline','Notch','off')

% axis([0 40 0.004 30]); % for filtered only data
axis([0 10 0.0005 10]); % for filtered only data
set(gca,'YScale','log');
set(gca,'ytick',[0.001,0.01,0.1,1])
set(gca,'xtick',1.55: 2.28 : 10)
set(gca,'xticklabel',box_name_cell)
xtickangle(90);

% xlabel('Time (s)','FontSize',12,'Color','k'); 
ylabel('Mean Z score','FontSize',12,'Color','k'); 

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 7, 'FontName','Arial'); %grid on;

    if section_num == 1; % title('\fontsize{7}Activity of Baseline cells');
elseif section_num == 2; % title('\fontsize{7}Activity of CS cells');
elseif section_num == 3; % title('\fontsize{7}Activity of US cells');
elseif section_num == 4; % title('\fontsize{7}Activity of ITI-1 cells');
% elseif section_num == 5; title('\fontsize{7}Activity of ITI-2 cells');
% % elseif section_num == 6; title('\fontsize{7}Activity of Rest cells');
% % elseif section_num == 7; title('\fontsize{7}Activity of HC-2 cells');
% % elseif section_num == 8; title('\fontsize{7}Activity of HC-3 cells');
% % elseif section_num == 7; title('\fontsize{7}Activity of STM(Acc) cells)');
% % elseif section_num == 8; title('\fontsize{7}Activity of STM(CS) cells');
% % elseif section_num == 9; title('\fontsize{7}Activity of STM(ITI) cells');
% elseif section_num == 6; title('\fontsize{7}Activity of LTM(Acc) cells');
% elseif section_num == 7; % title('\fontsize{7}Activity of LTM(CS) cells');
% elseif section_num == 8; title('\fontsize{7}Activity of LTM(ITI) cells');

    end
%%

elseif input_session_num == 9
% section_num = 2;
Box_X = [];
Box_Y = [];


for i = 1:input_session_num
        cat_x_temp_G1 = (cell2mat(G1_results_ratio(i,section_num)))';
        cat_x_temp_G2 = (cell2mat(G2_results_ratio(i,section_num)))';
    Box_X = [Box_X; cat_x_temp_G1];
    Box_X = [Box_X; cat_x_temp_G2];  
    
             if i == 1; group1 = ('Baseline');  group2 = ('Baseline-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 2; group1 = ('CS');  group2 = ('CS-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 3; group1 = ('US');  group2 = ('US-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];   
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 4; group1 = ('ITI-1');  group2 = ('ITI-1-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 5; group1 = ('ITI-2');  group2 = ('ITI-2-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 6; group1 = ('Rest-1');  group2 = ('Rest-1-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 7; group1 = ('HC-2');  group2 = ('HC-2-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 8; group1 = ('HC-3');  group2 = ('HC-3-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 7; group1 = ('STM-Baseline');  group2 = ('STM-Baseline-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 8; group1 = ('STM-CS');  group2 = ('STM-CS-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 9; group1 = ('STM-ITI');  group2 = ('STM-ITI-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 7; group1 = ('LTM-Baseline');  group2 = ('LTM-baseline-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 8; group1 = ('LTM-CS'); group2 = ('LTM-CS-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 9; group1 = ('LTM-ITI'); group2 = ('LTM-ITI-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];    
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         end  
    end

%% Boxplot new
% figure('color',[1,1,1]);
figure('Position',[244,181,180,130]); %[left bottom width height]
boxplot(Box_X, Box_Y,'factorgap',2,'color','kr','OutlierSize',2, 'Widths' ,1.0,'Symbol','','BoxStyle','outline','Notch','off')

% axis([0 40 0.004 30]); % for filtered only data
axis([0 25 -0.05 0.65]); % for filtered only data
% set(gca,'YScale','log');
set(gca,'xtick',1.55: 2.7 : 25)
set(gca,'xticklabel',{'Baseline', 'CS','US','ITI-1','ITI-2','Rest','LTM(Acc)','LTM(CS)','LTM(ITI)'})
xtickangle(90);

% xlabel('Time (s)','FontSize',12,'Color','k'); 
ylabel('Mean Z score','FontSize',12,'Color','k'); 

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 7, 'FontName','Arial'); %grid on;

    if section_num == 1; % title('\fontsize{7}Activity of Baseline cells');
elseif section_num == 2; % title('\fontsize{7}Activity of CS cells');
elseif section_num == 3; % title('\fontsize{7}Activity of US cells');
elseif section_num == 4; title('\fontsize{7}Activity of ITI-1 cells');
elseif section_num == 5; title('\fontsize{7}Activity of ITI-2 cells');
elseif section_num == 6; title('\fontsize{7}Activity of Rest cells');
% elseif section_num == 7; title('\fontsize{7}Activity of HC-2 cells');
% elseif section_num == 8; title('\fontsize{7}Activity of HC-3 cells');
% elseif section_num == 7; title('\fontsize{7}Activity of STM(Acc) cells)');
% elseif section_num == 8; title('\fontsize{7}Activity of STM(CS) cells');
% elseif section_num == 9; title('\fontsize{7}Activity of STM(ITI) cells');
elseif section_num == 7; title('\fontsize{7}Activity of LTM(Acc) cells');
elseif section_num == 8; % title('\fontsize{7}Activity of LTM(CS) cells');
elseif section_num == 9; title('\fontsize{7}Activity of LTM(ITI) cells');
    end
%%

%%
elseif input_session_num == 12
% section_num = 2;
Box_X = [];
Box_Y = [];


for i = 2:input_session_num
        cat_x_temp_G1 = (cell2mat(G1_results_ratio(i,section_num)))';
        cat_x_temp_G2 = (cell2mat(G2_results_ratio(i,section_num)))';
    Box_X = [Box_X; cat_x_temp_G1];
    Box_X = [Box_X; cat_x_temp_G2];  
    
             if i == 2; group1 = ('CS');  group2 = ('CS-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 3; group1 = ('US');  group2 = ('US-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];   
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 4; group1 = ('ITI-1');  group2 = ('ITI-1-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 5; group1 = ('ITI-2');  group2 = ('ITI-2-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 6; group1 = ('Rest-1');  group2 = ('Rest-1-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 7; group1 = ('HC-2');  group2 = ('HC-2-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 8; group1 = ('HC-3');  group2 = ('HC-3-g2');
%      Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 7; group1 = ('STM-Baseline');  group2 = ('STM-Baseline-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 8; group1 = ('STM-CS');  group2 = ('STM-CS-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 9; group1 = ('STM-ITI');  group2 = ('STM-ITI-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 10; group1 = ('LTM-Baseline');  group2 = ('LTM-baseline-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 11; group1 = ('LTM-CS'); group2 = ('LTM-CS-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 12; group1 = ('LTM-ITI'); group2 = ('LTM-ITI-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];    
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         end  
    end

%% Boxplot new
figure('color',[1,1,1]);
boxplot(Box_X, Box_Y,'factorgap',2,'color','kr','OutlierSize',2, ...
    'Widths' ,1.0,'Symbol','','BoxStyle','outline','Notch','off')

% axis([0 40 0.004 30]); % for filtered only data
axis([0 32 -0.1 0.6]); % for filtered only data
% set(gca,'YScale','log');
set(gca,'xtick',1.55: 2.85 : 32)
set(gca,'xticklabel',{'CS','US','ITI-1','ITI-2','Rest','STM(Acc)','STM(CS)','STM(ITI)',...
        'LTM(Acc)','LTM(CS)','LTM(ITI)'})
xtickangle(90);

xlabel('Time (s)','FontSize',12,'FontWeight','bold','Color','k'); 
ylabel('Z-score mean','FontSize',12,'FontWeight','bold','Color','k'); 

    if section_num == 2; title('\fontsize{7}Activity of CS cells');
elseif section_num == 3; title('\fontsize{7}Activity of US cells');
elseif section_num == 4; title('\fontsize{7}Activity of ITI-1 cells');
elseif section_num == 5; title('\fontsize{7}Activity of ITI-2 cells');
elseif section_num == 6; title('\fontsize{7}Activity of Rest cells');
% elseif section_num == 7; title('\fontsize{7}Activity of HC-2 cells');
% elseif section_num == 8; title('\fontsize{7}Activity of HC-3 cells');
elseif section_num == 7; title('\fontsize{7}Activity of STM(Acc) cells)');
elseif section_num == 8; title('\fontsize{7}Activity of STM(CS) cells');
elseif section_num == 9; title('\fontsize{7}Activity of STM(ITI) cells');
elseif section_num == 10; title('\fontsize{7}Activity of LTM(Acc) cells');
elseif section_num == 11; title('\fontsize{7}Activity of LTM(CS) cells');
elseif section_num == 12; title('\fontsize{7}Activity of LTM(ITI) cells');
    end
%%

elseif input_session_num == 14
% section_num = 2;
Box_X = [];
Box_Y = [];


for i = 2:input_session_num
        cat_x_temp_G1 = (cell2mat(G1_results_ratio(i,section_num)))';
        cat_x_temp_G2 = (cell2mat(G2_results_ratio(i,section_num)))';
    Box_X = [Box_X; cat_x_temp_G1];
    Box_X = [Box_X; cat_x_temp_G2];  
    
             if i == 2; group1 = ('CS');  group2 = ('CS-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 3; group1 = ('US');  group2 = ('US-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];   
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 4; group1 = ('ITI-1');  group2 = ('ITI-1-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 5; group1 = ('ITI-2');  group2 = ('ITI-2-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 6; group1 = ('HC-1');  group2 = ('HC-1-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 7; group1 = ('HC-2');  group2 = ('HC-2-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 8; group1 = ('HC-3');  group2 = ('HC-3-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 9; group1 = ('STM-Baseline');  group2 = ('STM-Baseline-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 10; group1 = ('STM-CS');  group2 = ('STM-CS-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 11; group1 = ('STM-ITI');  group2 = ('STM-ITI-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 12; group1 = ('LTM-Baseline');  group2 = ('LTM-baseline-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 13; group1 = ('LTM-CS'); group2 = ('LTM-CS-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];     
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         elseif i == 14; group1 = ('LTM-ITI'); group2 = ('LTM-ITI-g2');
     Box_Y_temp = (repmat({group1},[size(cell2mat(G1_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];    
     Box_Y_temp = (repmat({group2},[size(cell2mat(G2_results_ratio(i,section_num))),1]))'; Box_Y = [Box_Y; Box_Y_temp];
         end  
    end

%% Boxplot new
figure('color',[1,1,1]);
boxplot(Box_X, Box_Y,'factorgap',2,'color','kr','OutlierSize',2, ...
    'Widths' ,1.0,'Symbol','','BoxStyle','outline','Notch','off')

% axis([0 40 0.004 30]); % for filtered only data
axis([0 40 0.004 30]); % for filtered only data
% set(gca,'YScale','log');
set(gca,'xtick',1.85: 3 : 39)
set(gca,'xticklabel',{'CS','US','ITI-1','ITI-2','HC-1','HC-2','HC-3','STM(Acc)','STM(CS)','STM(ITI)',...
        'LTM(Acc)','LTM(CS)','LTM(ITI)'})
xtickangle(90);

xlabel('Time (s)','FontSize',12,'FontWeight','bold','Color','k'); 
ylabel('Responsiveness (vs. Baseline)','FontSize',12,'FontWeight','bold','Color','k'); 

    if section_num == 2; title('\fontsize{7}Activity of CS cells');
elseif section_num == 3; title('\fontsize{7}Activity of US cells');
elseif section_num == 4; title('\fontsize{7}Activity of ITI-1 cells');
elseif section_num == 5; title('\fontsize{7}Activity of ITI-2 cells');
elseif section_num == 6; title('\fontsize{7}Activity of HC-1 cells');
elseif section_num == 7; title('\fontsize{7}Activity of HC-2 cells');
elseif section_num == 8; title('\fontsize{7}Activity of HC-3 cells');
elseif section_num == 9; title('\fontsize{7}Activity of STM(Acc) cells)');
elseif section_num == 10; title('\fontsize{7}Activity of STM(CS) cells');
elseif section_num == 11; title('\fontsize{7}Activity of STM(ITI) cells');
elseif section_num == 12; title('\fontsize{7}Activity of LTM(Acc) cells');
elseif section_num == 13; title('\fontsize{7}Activity of LTM(CS) cells');
elseif section_num == 14; title('\fontsize{7}Activity of LTM(ITI) cells');
    end
%%
end

end