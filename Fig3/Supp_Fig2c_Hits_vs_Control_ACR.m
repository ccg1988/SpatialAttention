% by CCG @ 2025-07-29

clear; clc; 
% close all
load('SRF_shift_test_analyses_control.mat')

HvP_CL = SRF_shift_test_analyses_control(2).hits_vs_control; % *
increase_sig_CL = SRF_shift_test_analyses_control(2).sig_increasesHvC_index;
decrease_sig_CL = SRF_shift_test_analyses_control(2).sig_decreasesHvC_index;
HvP_A1 = SRF_shift_test_analyses_control(1).hits_vs_control; % square
increase_sig_A1 = SRF_shift_test_analyses_control(1).sig_increasesHvC_index;
decrease_sig_A1 = SRF_shift_test_analyses_control(1).sig_decreasesHvC_index;
HvP_R = SRF_shift_test_analyses_control(3).hits_vs_control; % >
increase_sig_R = SRF_shift_test_analyses_control(3).sig_increasesHvC_index;
decrease_sig_R = SRF_shift_test_analyses_control(3).sig_decreasesHvC_index;

pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position',[X_size*0.1 Y_size*0.3 X_size*0.3 Y_size*0.3]);
sz = 20 ;

scatter(HvP_CL(increase_sig_CL,1), HvP_CL(increase_sig_CL,2), sz*1.5, '*', 'MarkerEdgeColor', rgb('BlueViolet'),'LineWidth', 1); hold on
scatter(HvP_CL(decrease_sig_CL,1), HvP_CL(decrease_sig_CL,2), sz*1.5, '*', 'MarkerEdgeColor', rgb('Green'),'LineWidth', 1);
scatter(HvP_A1(increase_sig_A1,1), HvP_A1(increase_sig_A1,2), sz*1.5, 's', 'MarkerEdgeColor', rgb('BlueViolet'),'LineWidth', 1); 
scatter(HvP_A1(decrease_sig_A1,1), HvP_A1(decrease_sig_A1,2), sz*1.5, 's', 'MarkerEdgeColor', rgb('Green'),'LineWidth', 1);
scatter(HvP_R(increase_sig_R,1), HvP_R(increase_sig_R,2), sz*1.5, '>', 'MarkerEdgeColor', rgb('BlueViolet'),'LineWidth', 1); 
scatter(HvP_R(decrease_sig_R,1), HvP_R(decrease_sig_R,2), sz*1.5, '>', 'MarkerEdgeColor', rgb('Green'),'LineWidth', 1);
plot(-15:105, -15:105, 'k--');
xlim([-15, 105]); ylim([-15, 105]); 
title({['A1 total=',num2str(length(increase_sig_A1)), ' inc=', num2str(length(find(increase_sig_A1))), ' dec=', num2str(length(find(decrease_sig_A1)))], ...
    ['CL total=',num2str(length(increase_sig_CL)), ' inc=', num2str(length(find(increase_sig_CL))), ' dec=', num2str(length(find(decrease_sig_CL)))], ...
     ['R total=',num2str(length(increase_sig_R)), ' inc=', num2str(length(find(increase_sig_R))), ' dec=', num2str(length(find(decrease_sig_R)))]})
pbaspect([1 1 1])
disp('A1 CL R: Number and percent of significant changes (relative to all locations)')
disp([(66+7)/247, (38+5)/120, (11+11)/120])
disp('A1 CL R: Number and percent of significant increase (relative to significant locations)')
disp([66/(66+7), 38/(38+5), 11/(11+11)])
