function fxnHF_1s_bin_matrix(ca_1s_bin1, ca_1s_bin2, ca_1s_bin_mean1, ca_1s_bin_mean2, reshape_range)
% rehsape_range = (151:600)

 %% ttest load data
 %[h,p] = ttest(tx,ty)

 for i = 1:size(ca_1s_bin1,1)  
    tx = ca_1s_bin1(i,:);  
    ty = ca_1s_bin2(i,:); 
    [ttest_table(i,1),ttest_table(i,2)] = ttest2(tx,ty);
 end
    %% reshape
    
%     rehsape_range = (151:600);
    
    rs_temp1 = ca_1s_bin_mean1(reshape_range);
    rs_temp2 = ca_1s_bin_mean2(reshape_range);
    rs_temp3 = ttest_table(reshape_range,2);
        
        rs_mean1 = reshape(rs_temp1,45,[]);
        rs_mean2 = reshape(rs_temp2,45,[]);
        rs_mean3 = reshape(rs_temp3,45,[]);
        
        rs_mean1_2 = rs_mean1-rs_mean2;
        
        %ax1 = subplot(2,1,1); colormap(ax1,spring)
        
%         figure; 
%         ax1 = subplot(411); imagesc(rs_mean1')   ; colormap(ax1,parula);  clim([0 0.2]);         
%             c = colorbar;  c.Label.String = 'Mean z-score'; xlabel('Time (s)'); ylabel('Trial#');   
%         ax2 = subplot(412); imagesc(rs_mean2')   ; colormap(ax2,parula);  clim([0 0.2]);
%             c = colorbar;  c.Label.String = 'Mean z-score'; xlabel('Time (s)'); ylabel('Trial#');
%         ax3 = subplot(413); imagesc(rs_mean1_2') ; colormap(ax3,jet); clim([-0.25 0.25]);      
%             c = colorbar;  c.Label.String = 'Subtaction value (Cntr-KO)'; xlabel('Time (s)'); ylabel('Trial#');
%         ax4 = subplot(414); imagesc(rs_mean3')   ; clim([0.049999 0.5]);  colormap(ax4,summer);
%             c = colorbar;  c.Label.String = 't-test p-value'; xlabel('Time (s)'); ylabel('Trial#');
    %% reshape --> mean
    
%     rs_mean1_mean_15 = mean(rs_mean1(:,1:5),2);
%     rs_mean2_mean_15 = mean(rs_mean2(:,1:5),2);
%     rs_mean3_mean_15 = mean(rs_mean3(:,1:5),2);
%     rs_mean1_2_mean_15 = mean(rs_mean1_2(:,1:5),2);
    
%     figure;
%     plot(rs_mean1_mean_15);
%     hold on
%     plot(rs_mean2_mean_15);
%     hold on
%     plot(rs_mean3_mean_15);
%     hold on
%     plot(rs_mean1_2_mean_15);
    
%     %%
%     rs_mean1_mean_610 = mean(rs_mean1(:,6:10),2);
%     rs_mean2_mean_610 = mean(rs_mean2(:,6:10),2);
%     rs_mean3_mean_610 = mean(rs_mean3(:,6:10),2);
%     rs_mean1_2_mean_610 = mean(rs_mean1_2(:,6:10),2);
    
%     figure;
%     plot(rs_mean1_mean_610);
%     hold on
%     plot(rs_mean2_mean_610);
%     hold on
%     plot(rs_mean3_mean_610);
%     hold on
%     plot(rs_mean1_2_mean_610);
%     
%     
%     legend('1','2','3','4');
 

%% all trial normarization
% clear bin_temp_ii bin_temp_rs_mean1 bin_temp_rs_mean2
plot_range = reshape_range(136:end);

    bin_temp1 = ca_1s_bin1;
    bin_temp2 = ca_1s_bin2;
    
    for ii = 1:size(bin_temp1,2) 
        bin_temp_ii = bin_temp1(plot_range,ii);    
        bin_temp_rs_mean1(:,ii) = mean(reshape(bin_temp_ii,45,[]),2);
    end
        bin_temp_rs_mean1 = bin_temp_rs_mean1';
        
    for ii = 1:size(bin_temp2,2) 
        bin_temp_ii = bin_temp2(plot_range,ii);    
        bin_temp_rs_mean2(:,ii) = mean(reshape(bin_temp_ii,45,[]),2);
    end
        bin_temp_rs_mean2 = bin_temp_rs_mean2';
        
        bin_all_table(1,:) = mean(bin_temp_rs_mean1(:,:),1);
        bin_all_table(4,:) = std(bin_temp_rs_mean1(:,:),0,1) ./ sqrt(size(bin_temp_rs_mean1,1));
        
        bin_all_table(2,:) = mean(bin_temp_rs_mean2(:,:),1);
        bin_all_table(5,:) = std(bin_temp_rs_mean2(:,:),0,1) ./ sqrt(size(bin_temp_rs_mean2,1));
     
% ttest load data
% [h,p] = ttest(tx,ty)
    
    for iii = 1:size(bin_all_table,2)
    
    bin_tx = bin_temp_rs_mean1(:,iii);  
    bin_ty = bin_temp_rs_mean2(:,iii); 
    
    [bin_all_table(7,iii), bin_all_table(8,iii)] = ttest2(bin_tx,bin_ty);
    
    end
% figure, % errorbar(y, error)
%        
%     figure;
%     errorbar(bin_all_table(1,:), bin_all_table(4,:));
%     hold on
%     errorbar(bin_all_table(2,:), bin_all_table(5,:));
%     title('Trial ALL normalized');
%%
figure('Position',[244,181,200,700]); %[left bottom width height]
ax1 = subplot(511); imagesc(rs_mean1')   ; colormap(ax1,parula);  clim([0 0.35]);         
c.Label.String = 'Mean z-score'; xlabel('Time (s)'); % ylabel('Trial#'); 
set(gca,'ytick',1:13); set(gca,'yticklabel',{'(Pre -135s)',' ','(Pre -45s)',' ','2',' ','4',' ', '6', ' ','8',' ','10',})
title('1 s-binned z-score (Ctrl)');

ax2 = subplot(512); imagesc(rs_mean2')   ; colormap(ax2,parula);  clim([0 0.35]);
c.Label.String = 'Mean z-score'; xlabel('Time (s)'); % ylabel('Trial#');
set(gca,'ytick',1:13); set(gca,'yticklabel',{'(Pre -135s)',' ','(Pre -45s)',' ','2',' ','4',' ', '6', ' ','8',' ','10',})
title('1 s-binned z-score (KO)');

subplot(5,1,3); %'--or'
errorbar(bin_all_table(1,:), bin_all_table(4,:),'-ok','MarkerSize',3,'CapSize',0,'MarkerEdgeColor','k','MarkerFaceColor','w'); 
hold on
errorbar(bin_all_table(2,:), bin_all_table(5,:),'-or','MarkerSize',3,'CapSize',0,'MarkerEdgeColor','r','MarkerFaceColor','w');
xlim([0 45]); title('Normalized z-score over 10 trials'); xlabel('Time (s)'); ylabel('Mean z-score'); grid on; %legend('Ctrl','KO')
set(gca,'YScale','log'); set(gca,'ytick',[0.01, 0.1, 1]);
ylim([0.01 1.0]);

% for colorbar picking up
ax1 = subplot(514); imagesc(rs_mean1')   ; colormap(ax1,parula);  clim([0 0.35]);         
c = colorbar;  c.Label.String = 'Mean z-score'; xlabel('Time (s)'); % ylabel('Trial#'); 
set(gca,'ytick',1:13)
set(gca,'yticklabel',{'(Pre -135s)',' ','(Pre -45s)',' ','2',' ','4',' ', '6', ' ','8',' ','10',})

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 7, 'FontName','Arial')
            %%
end