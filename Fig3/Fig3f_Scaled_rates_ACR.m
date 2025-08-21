clear; clc; 
close all
load('SRF_shift_test_analyses.mat')
makers_all = {'s','*','>'}; % A1 CL R
pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position',[X_size*0.1 Y_size*0.3 X_size*0.3 Y_size*0.3]);
sz = 45 ;
HmP_CL = [];
HmP_A1R = [];
for p = 1 : 3
hits_vs_passive = SRF_shift_test_analyses(p).hits_vs_passive;

scaled_value_1 = hits_vs_passive(:, data_columns.scaled_value_1);
scaled_value_2 = hits_vs_passive(:, data_columns.scaled_value_2); % value_3 & scaled_value_3 are always empty

sig_increase = hits_vs_passive(:, data_columns.sig_increase); % sig_difference sig_increase sig_decrease
sig_increase = logical(sig_increase); 

P = scaled_value_1(sig_increase);
H = scaled_value_2(sig_increase);
HmP = H-P;
if p == 2
    HmP_CL = HmP;
else
    HmP_A1R = [HmP_A1R; HmP];
end
scatter(P(HmP>1), HmP(HmP>1), sz*1, makers_all{p}, 'MarkerEdgeColor', rgb('Purple'),'LineWidth', 1); hold on
scatter(P(HmP<1), HmP(HmP<1), sz*1, makers_all{p}, 'MarkerEdgeColor', rgb('Violet'),'LineWidth', 1);
scatter(median(P), median(HmP), sz*1.5, makers_all{p}, 'MarkerEdgeColor', rgb('Black'),'LineWidth', 2)
yline(1, 'k--'); hold on
xlim([0, 1]); xticks([0 0.2 0.4 0.6 0.8 1]); 
ylim([0, 4]); yticks([0 1 2 3 4]);
xlabel('Scaled rate passive'); ylabel('Scaled increase hits')
per = num2str(round(100*(length(find(HmP>1))/length(P))));
disp(['Hits minus Passive median=',num2str(median(HmP)), ' Passive median=',num2str(median(P))])
disp(['Total=',num2str(length(P)), ' >1=',num2str(length(find(HmP>1))), ' percent=', per, '%'])
pbaspect([1 1 1])
end
disp(['Hits minus Passive CL vs A1&R=', num2str(ranksum(HmP_CL, HmP_A1R))])