% by CCG @ 2025-07-29
% Select all the 487 target locations (158 units) from Fig3a that was significantly driven by sound under either passive or hits
% However, the total locations

clear; clc; 
% close all
load('SRF_shift_test_analyses.mat')

HvP_CL = SRF_shift_test_analyses(2).hits_vs_passive;
increase_sig_CL = HvP_CL(:, data_columns.sig_increase)==1;
decrease_sig_CL = HvP_CL(:, data_columns.sig_decrease)==1;
HvP_A1 = SRF_shift_test_analyses(1).hits_vs_passive;
increase_sig_A1 = HvP_A1(:, data_columns.sig_increase)==1;
decrease_sig_A1 = HvP_A1(:, data_columns.sig_decrease)==1;
HvP_R = SRF_shift_test_analyses(3).hits_vs_passive;
increase_sig_R = HvP_R(:, data_columns.sig_increase)==1;
decrease_sig_R = HvP_R(:, data_columns.sig_decrease)==1;

pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position',[X_size*0.1 Y_size*0.3 X_size*0.3 Y_size*0.3]);
sz = 50 ;

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

% Sound driven: 74+52+32 = 158(93+65) units
% Hits modulation: significant Inc=70, significant Dec=19 
HvP_units_A1 = SRF_shift_test_analyses(1).unitsHvP; % 74=42+32  Sig Inc&Dec=36&6
HvP_units_CL = SRF_shift_test_analyses(2).unitsHvP; % 52=38+14  Sig Inc&Dec=24&5
HvP_units_R = SRF_shift_test_analyses(3).unitsHvP;   % 32=13+19  Sig Inc&Dec=10&8

% figure;
% plot(HvP_CL(:,data_columns.is_driven1)); hold on
% plot(HvP_CL(:,data_columns.is_driven2))
% plot(HvP_CL(:,data_columns.is_driven1)+HvP_CL(:,data_columns.is_driven2), 'LineWidth', 2) % columns 13=Hits and 14=Passive, 15=Control is empty
% legend({'Hits', 'Passive', 'Added'})
% title('All the location/row was either driven by Hits or Passive')