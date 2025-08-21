clear; clc; 
close all
CH = 4 ;
stim_select = 3 ;             % relative stimulus in stim_list; <=5
stim_select_P = 7 ;         % relative=fixed stimulus in stim_list
bin_ms = 50 ;

file_name='M9X0832';
[spikeTimes, spikeMat, pre_ms, stim_ms, post_ms, t_ms, stim_list] = spike_extract(file_name, CH);
nrep = numel(spikeTimes{stim_select});
file_name_P='M9X0842';
[spikeTimes_P, spikeMat_P, ~, ~, ~, ~, stim_list_P] = spike_extract(file_name_P, CH);

pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position',[X_size*0.05 Y_size*0.2 X_size*0.6 Y_size*0.3]);

subplot(1,3,1)
patch([0 stim_ms stim_ms 0], [0 0 nrep+1 nrep+1], [0.9 0.9 0.9],'EdgeColor', 'none', 'FaceAlpha', 0.25); hold on
for r = 1:nrep
    x = spikeTimes{stim_select}{r};
    if ~isempty(x)
        plot(x, r*ones(size(x)), 'k.', 'MarkerSize', 8); hold on
    end
end
xlim([t_ms(1) t_ms(end)]); ylim([0 nrep + 1]); xlabel('Time from stimulus onset (ms)'); ylabel('Repeat #');
title(sprintf('Hits Location %d', stim_list(stim_select))); set(gca, 'Box', 'off');

subplot(1,3,2)
nrep = numel(spikeTimes_P{stim_select_P});
patch([0 stim_ms stim_ms 0], [0 0 nrep+1 nrep+1], [0.9 0.9 0.9],'EdgeColor', 'none', 'FaceAlpha', 0.25);hold on
for r = 1:nrep
    x = spikeTimes_P{stim_select_P}{r};
    if ~isempty(x)
        plot(x, r*ones(size(x)), 'k.', 'MarkerSize', 8); hold on
    end
end
xlim([t_ms(1) t_ms(end)]); ylim([0 nrep + 1]); xlabel('Time from stimulus onset (ms)'); ylabel('Repeat #');
title(sprintf('Passive Location %d', stim_list_P(stim_select_P))); set(gca, 'Box', 'off');
%%
[t_msX, mean_rateX_bs, baseline_rate] = psth_plot(spikeMat, stim_select, bin_ms, pre_ms, stim_ms, post_ms, t_ms);
[~, mean_rateX_bs_P, baseline_rate_P] = psth_plot(spikeMat_P, stim_select_P, bin_ms, pre_ms, stim_ms, post_ms, t_ms);

% pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
% figure('position',[X_size*0.05 Y_size*0.2 X_size*0.45 Y_size*0.3]);
subplot(1,3,3)

yMax = max(mean_rateX_bs);
if yMax <= 0, yMax = eps; end
patch([0 stim_ms stim_ms 0], [0 0 yMax yMax], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.25); hold on
plot(t_msX, mean_rateX_bs, 'r', 'LineWidth', 1.5); hold on % stairs
plot(t_msX, mean_rateX_bs_P, 'b', 'LineWidth', 1.5);

yline(0, '--', 'Color', 'k');
legend({'', 'Hits', 'Passive'})

xlabel('Time from stimulus onset (ms)'); ylabel('Firing rate (Hz, minus-spont)');
title(['Spont. rate Hits=', num2str(round(baseline_rate,2)), ' Passive=', num2str(round(baseline_rate_P,2))])
box off; axis tight;


