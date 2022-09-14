function fxn_boxplot_double_for_ca3_revised_importdata(input_cell, yaxis_max, comparison_section, ctrl_range, ko_range);
%%
    % section_num = 2;
Box_X = [];
Box_Y = [];
Box_indiv_sort = [];
   
%%
if comparison_section == 1
for i = 1:9
cat_x_temp_G1 = ((input_cell(ctrl_range,i)));
cat_x_temp_G2 = ((input_cell(ko_range,i)));
Box_X = [Box_X; cat_x_temp_G1];
Box_X = [Box_X; cat_x_temp_G2];  


if i == 1; group1 = ('Baseline');  group2 = ('Baseline-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 2; group1 = ('CS');  group2 = ('CS-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 3; group1 = ('US');  group2 = ('US-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];   
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 4; group1 = ('ITI-1');  group2 = ('ITI-1-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 5; group1 = ('ITI-2');  group2 = ('ITI-2-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 6; group1 = ('Rest-1');  group2 = ('Rest-1-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
%          elseif i == 7; group1 = ('HC-2');  group2 = ('HC-2-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 8; group1 = ('HC-3');  group2 = ('HC-3-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 7; group1 = ('STM-Baseline');  group2 = ('STM-Baseline-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 8; group1 = ('STM-CS');  group2 = ('STM-CS-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 9; group1 = ('STM-ITI');  group2 = ('STM-ITI-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
elseif i == 7; group1 = ('LTM-Baseline');  group2 = ('LTM-baseline-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 8; group1 = ('LTM-CS'); group2 = ('LTM-CS-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 9; group1 = ('LTM-ITI'); group2 = ('LTM-ITI-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];    
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
end
end
%%
elseif comparison_section == 2

for i = 1:9
cat_x_temp_G1 = ((input_cell(ctrl_range,11+i)));
cat_x_temp_G2 = ((input_cell(ko_range,11+i)));
Box_X = [Box_X; cat_x_temp_G1];
Box_X = [Box_X; cat_x_temp_G2];  

if i == 1; group1 = ('Baseline');  group2 = ('Baseline-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 2; group1 = ('CS');  group2 = ('CS-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 3; group1 = ('US');  group2 = ('US-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];   
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 4; group1 = ('ITI-1');  group2 = ('ITI-1-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 5; group1 = ('ITI-2');  group2 = ('ITI-2-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 6; group1 = ('Rest-1');  group2 = ('Rest-1-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
%          elseif i == 7; group1 = ('HC-2');  group2 = ('HC-2-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 8; group1 = ('HC-3');  group2 = ('HC-3-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 7; group1 = ('STM-Baseline');  group2 = ('STM-Baseline-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 8; group1 = ('STM-CS');  group2 = ('STM-CS-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 9; group1 = ('STM-ITI');  group2 = ('STM-ITI-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
elseif i == 7; group1 = ('LTM-Baseline');  group2 = ('LTM-baseline-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 8; group1 = ('LTM-CS'); group2 = ('LTM-CS-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 9; group1 = ('LTM-ITI'); group2 = ('LTM-ITI-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];    
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
end
end
%%
elseif comparison_section == 3

for i = 1:9
cat_x_temp_G1 = ((input_cell(ctrl_range,23+i)));
cat_x_temp_G2 = ((input_cell(ko_range,23+i)));
Box_X = [Box_X; cat_x_temp_G1];
Box_X = [Box_X; cat_x_temp_G2];  

if i == 1; group1 = ('Baseline');  group2 = ('Baseline-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 2; group1 = ('CS');  group2 = ('CS-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 3; group1 = ('US');  group2 = ('US-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];   
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 4; group1 = ('ITI-1');  group2 = ('ITI-1-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 5; group1 = ('ITI-2');  group2 = ('ITI-2-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 6; group1 = ('Rest-1');  group2 = ('Rest-1-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
%          elseif i == 7; group1 = ('HC-2');  group2 = ('HC-2-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 8; group1 = ('HC-3');  group2 = ('HC-3-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 7; group1 = ('STM-Baseline');  group2 = ('STM-Baseline-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 8; group1 = ('STM-CS');  group2 = ('STM-CS-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
%          elseif i == 9; group1 = ('STM-ITI');  group2 = ('STM-ITI-g2');
%      Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
%      Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
elseif i == 7; group1 = ('LTM-Baseline');  group2 = ('LTM-baseline-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 8; group1 = ('LTM-CS'); group2 = ('LTM-CS-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];     
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
elseif i == 9; group1 = ('LTM-ITI'); group2 = ('LTM-ITI-g2');
Box_Y_temp = (repmat({group1},[size((input_cell(ctrl_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];    
Box_Y_temp = (repmat({group2},[size((input_cell(ko_range,11+i))),1])); Box_Y = [Box_Y; Box_Y_temp];
indiv_tmep = (repmat( ((i*2)-1) ,[size((input_cell(ctrl_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
indiv_tmep = (repmat( (i*2)     ,[size((input_cell(ko_range,i))),1])); Box_indiv_sort = [Box_indiv_sort; indiv_tmep];
end
end
end
%% Boxplot new
% figure('color',[1,1,1]);
figure('Position',[244,181,172,130]); %[left bottom width height]
boxplot(Box_X, Box_Y,'factorgap',2,'color','kr','OutlierSize',2, 'Widths' ,1.0,'Symbol','','BoxStyle','outline','Notch','off')
% boxplot(Box_X, Box_Y,'factorgap',2,'color','kr','OutlierSize',2, 'Widths' ,1.0,'Symbol','','BoxStyle','outline','Notch','off') % outlier on

% 211129 add, individual plot
hold on;
scatter( (Box_indiv_sort(:).*1.35)-0.5 , Box_X(:), 1.5,'filled','MarkerFaceAlpha',0.7','jitter','on','jitterAmount',0.35);
hold off;

% axis([0 40 0.004 30]); % for filtered only data

% axis([0 25 0.00001 yaxis_max]); % for ca3 triple shared, log
axis([0 25 0 yaxis_max]); % for ca3 triple not-shared, normal

% set(gca,'YScale','log'); % for ca3 triple shared, log
set(gca,'xtick',1.55: 2.7 : 25)
set(gca,'xticklabel',{'Baseline', 'CS','US','ITI-1','ITI-2','Rest','LTM(Acc)','LTM(CS)','LTM(ITI)'})
xtickangle(90);

% xlabel('Time (s)','FontSize',12,'Color','k'); 
ylabel('Connectivity','FontSize',12,'Color','k'); 

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 7, 'FontName','Arial'); %grid on;

%     if input_session_num == 1; % title('\fontsize{7}Activity of Baseline cells');
% elseif input_session_num == 2; % title('\fontsize{7}Activity of CS cells');
% elseif input_session_num == 3; % title('\fontsize{7}Activity of US cells');
% elseif input_session_num == 4; title('\fontsize{7}Activity of ITI-1 cells');
% elseif input_session_num == 5; title('\fontsize{7}Activity of ITI-2 cells');
% elseif input_session_num == 6; title('\fontsize{7}Activity of Rest cells');
% % elseif input_session_num == 7; title('\fontsize{7}Activity of HC-2 cells');
% % elseif input_session_num == 8; title('\fontsize{7}Activity of HC-3 cells');
% % elseif input_session_num == 7; title('\fontsize{7}Activity of STM(Acc) cells)');
% % elseif input_session_num == 8; title('\fontsize{7}Activity of STM(CS) cells');
% % elseif input_session_num == 9; title('\fontsize{7}Activity of STM(ITI) cells');
% elseif input_session_num == 7; title('\fontsize{7}Activity of LTM(Acc) cells');
% elseif input_session_num == 8; % title('\fontsize{7}Activity of LTM(CS) cells');
% elseif input_session_num == 9; title('\fontsize{7}Activity of LTM(ITI) cells');
%     end
%%
end
