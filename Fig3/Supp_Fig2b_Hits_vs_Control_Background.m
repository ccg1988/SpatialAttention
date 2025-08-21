% by CCG @ 2025-07-28
% no location has both significant Increase and Decrease
% there are 1216 locations in 304 sessions since each session has 4 locations as Target
% and 1 location as Background (Fig3b). Each location has mutliple repeats.

clear; clc; close all
load('SRF_shift_test_analysis_control.mat')
AvP_units_id = SRF_shift_test_analysis_control.unitsBHvC; %************* Key Changes
N_units = size(AvP_units_id, 1);
rP_trial_mean = nan(N_units, 1);
rA_trial_mean = nan(N_units, 1);
rP_trial_sig = [] ;
rP_trial_nonsig = [] ;
rA_trial_sig = [] ;
rA_trial_nonsig = [] ;
AvP_units_raw = SRF_shift_test_analysis_control.background_hits_vs_control; %************* Key Changes
N_trials = size(AvP_units_raw, 1);
AvP_units = [] ;
id_rate_P = 1 ;
id_rate_A = 2 ;
id_sig = 19 ; 
id_sig_In = 20 ;
id_sig_De = 21 ;
id_mon = 30 ;
id_unit = 31 ;
for n = 1 : N_units
    mon_id = AvP_units_id(n, 1); unit_id = AvP_units_id(n, 2);
    temp_trial_id = find ((AvP_units_raw(:, id_mon)==mon_id)&(AvP_units_raw(:, id_unit)==unit_id));
    temp_trials = AvP_units_raw(temp_trial_id, :) ;
    rP_trial_mean(n) = mean(temp_trials(:, id_rate_P));   
    rA_trial_mean(n) = mean(temp_trials(:, id_rate_A));
%     rP_trial = [rP_trial; temp_trials(:, id_rate_P)];
%     rA_trial = [rA_trial; temp_trials(:, id_rate_A)];
%     if (AvP_units_raw(n, id_sig_In)+AvP_units_raw(n, id_sig_De))==1
%         AvP_units = [AvP_units; AvP_units_raw(n, :)] ;
%     end    
end    
MI_all=[];
MI_In=[];
MI_De=[];
for n = 1 : N_trials
    p_temp = AvP_units_raw(n, id_rate_P);
    a_temp = AvP_units_raw(n, id_rate_A);
    if AvP_units_raw(n, id_sig)==1 %>1 spike and significant different from spontaneous
        rP_trial_sig = [rP_trial_sig; p_temp];
        rA_trial_sig = [rA_trial_sig; a_temp];
        MI_temp = (a_temp - p_temp)/(a_temp + p_temp) ;
        if abs(MI_temp)>1
            MI_temp = sign(MI_temp);
        end     
        MI_all = [MI_all; MI_temp] ; %=SRF_shift_test_analysis.MI_hits_vs_passive_sig      
        if AvP_units_raw(n, id_sig_In)==1
            MI_In = [MI_In; MI_temp] ; %=SRF_shift_test_analysis.MI_hits_vs_passive_sig_increases
        elseif AvP_units_raw(n, id_sig_De)==1
            MI_De = [MI_De; MI_temp] ;
        end    
    elseif AvP_units_raw(n, id_sig)==0
        rP_trial_nonsig = [rP_trial_nonsig; p_temp];
        rA_trial_nonsig = [rA_trial_nonsig; a_temp];
    end    
end    
%%
pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position',[X_size*0.1 Y_size*0.3 X_size*0.3 Y_size*0.3]);
sz=30;

disp_min = min([rP_trial_sig; rA_trial_sig; rP_trial_nonsig; rA_trial_nonsig]);
disp_max = max([rP_trial_sig; rA_trial_sig; rP_trial_nonsig; rA_trial_nonsig]);
scatter(rP_trial_nonsig, rA_trial_nonsig, sz/2, 'o', 'MarkerEdgeColor', rgb('Silver'), 'MarkerFaceColor',rgb('White')); hold on
scatter(rP_trial_sig(rP_trial_sig<rA_trial_sig), rA_trial_sig(rP_trial_sig<rA_trial_sig), sz*1, 'x', 'MarkerEdgeColor', rgb('BlueViolet'),'LineWidth', 1); 
scatter(rP_trial_sig(rP_trial_sig>rA_trial_sig), rA_trial_sig(rP_trial_sig>rA_trial_sig), sz*1, 'x', 'MarkerEdgeColor', rgb('Green'),'LineWidth', 1); 
plot(disp_min:disp_max, disp_min:disp_max, 'k--'); hold on
xlim([disp_min disp_max]); xlabel('Passive Spikes/S')
ylim([disp_min disp_max]); ylabel('Hits Spikes/S')
title({['Non-significant stimuli above: ', num2str(length(find(rA_trial_nonsig>rP_trial_nonsig))), ' below: ', num2str(length(find(rA_trial_nonsig<rP_trial_nonsig)))], ...
    ['Significant stimuli above: ', num2str(length(find(rA_trial_sig>rP_trial_sig))), ' below: ', num2str(length(find(rA_trial_sig<rP_trial_sig)))]})
pbaspect([1 1 1])
