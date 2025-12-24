% by CCG @ 2025-12-23

clear; clc; close all
load('SRF_shift_test_analysis.mat')
AvP_units_id = SRF_shift_test_analysis.unitsHvP;
N_units = size(AvP_units_id, 1); %>>>>>>>>>>>>>>208
rP_trial_mean = nan(N_units, 1);
rA_trial_mean = nan(N_units, 1);
AvsP_p_value = nan(N_units, 1);
AvP_units_raw = SRF_shift_test_analysis.hits_vs_passive;
id_rate_P = 1 ;
id_rate_A = 2 ;
id_mon = 30 ;
id_unit = 31 ;
for n = 1 : N_units
    mon_id = AvP_units_id(n, 1); unit_id = AvP_units_id(n, 2);
    temp_trial_id = find ((AvP_units_raw(:, id_mon)==mon_id)&(AvP_units_raw(:, id_unit)==unit_id));
    temp_trials = AvP_units_raw(temp_trial_id, :) ;
    rP_trial_mean(n) = mean(temp_trials(:, id_rate_P));   
    rA_trial_mean(n) = mean(temp_trials(:, id_rate_A));
    [~, AvsP_p_value(n),~,~] = ttest(temp_trials(:, 1), temp_trials(:, 2)); % only one neuron (#78) have p-value of NaN (only 1 trial)
end    

%%
pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position',[X_size*0.1 Y_size*0.3 X_size*0.3 Y_size*0.3]);
sz=30;

disp_min = -10;
disp_max = 70;
scatter(rP_trial_mean, rA_trial_mean, sz/2, 'o', 'MarkerEdgeColor', rgb('Silver'), 'MarkerFaceColor',rgb('White')); hold on
condition1 = (AvsP_p_value < 0.05) & (rP_trial_mean < rA_trial_mean);
scatter(rP_trial_mean(condition1), rA_trial_mean(condition1), sz*1, 'x', 'MarkerEdgeColor', rgb('BlueViolet'),'LineWidth', 1); 
condition2 = (AvsP_p_value < 0.05) & (rP_trial_mean > rA_trial_mean);
scatter(rP_trial_mean(condition2), rA_trial_mean(condition2), sz*1, 'x', 'MarkerEdgeColor', rgb('Green'),'LineWidth', 1); 
plot(disp_min:disp_max, disp_min:disp_max, 'k--'); hold on
xlim([disp_min disp_max]); xlabel('Passive Spikes/S')
ylim([disp_min disp_max]); ylabel('Hits Spikes/S')
xticks([-10 10 30 50 70])
xticklabels({'-10','10','30','50','70'})
yticks([-10 10 30 50 70])
yticklabels({'-10','10','30','50','70'})
title({
    ['Significant neuron above: ', num2str(length(find(condition1))), ' below: ', num2str(length(find(condition2)))], []})
pbaspect([1 1 1])