clear; clc; 
close all
CH = 4 ;
stim_select = 5 ;           % relative stimulus in stim_list; <=5
bin_ms = 10 ;
pre_ms    = 200;           % prestimulus window
stim_ms   = 200;           % stimulus duration
post_ms   = 300;           % poststimulus window
t_ms   = -pre_ms : (stim_ms + post_ms - 1);   % e.g., -200 ... 499  (700 bins); 1 ms resolution
edges  = [t_ms, t_ms(end) + 1];               % histcounts edges
stim_col  = 1;   % stimulus ID
label_col = 3;   % event code: 4=spike, 7=stimulus onset, 8/9=lick, 10=reward
time_col  = 4;   % time in microseconds

file_name='M9X0842';
data_raw=str2num(file_name);
data=data_raw.data;
S=data_raw.stimulus_ch1; 
time_ms = double(data(:, time_col)) / 1000; % Convert us to ms
is_spike = data(:, label_col) == CH ; % Event masks
is_stim  = data(:, label_col) == 7 ;
stim_list = unique(S(:,2), 'stable') ;   % five stimuli and first is always "background"
spikeTimes = cell(numel(stim_list), 1);  % cell{stim}{repeat} = vector of spike times (ms, relative to onset)
spikeMat   = cell(numel(stim_list), 1);  % cell{stim} = [nRepeats x 700] counts per ms
for si = 1:numel(stim_list)
    sid = stim_list(si);
    ids = S(S(:,2) == sid, 1);                      % map type (8 vs 3-7-15-17) -> its 11 true IDs
    onset_idx = find(is_stim & ismember(data(:, stim_col), ids));
    onset_t   = time_ms(onset_idx); % All onsets for this stimulus
    nrep      = numel(onset_t);

    spikeTimes{si} = cell(nrep, 1);
    spikeMat{si}   = zeros(nrep, numel(t_ms), 'double');

    for r = 1:nrep
        t0 = onset_t(r);  % onset time (ms)

        win_start = t0 - pre_ms;
        win_end   = t0 + stim_ms + post_ms;

        % Spikes inside the window
        spk_t   = time_ms(is_spike & time_ms >= win_start & time_ms < win_end);
        rel_spk = spk_t - t0;   % relative to onset (ms)

        spikeTimes{si}{r}   = rel_spk(:);
        spikeMat{si}(r, :)  = histcounts(rel_spk, edges);
    end
end

%% -------------------- VISUALIZE single stimulus --------------------

nrep = numel(spikeTimes{stim_select});

pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position',[X_size*0.05 Y_size*0.2 X_size*0.6 Y_size*0.3]); 

subplot(1,2, 1)
patch([0 stim_ms stim_ms 0], [0 0 nrep+1 nrep+1], [0.9 0.9 0.9], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.25); hold on
% Raster: each repeat is a row, spikes are points
for r = 1:nrep
    x = spikeTimes{stim_select}{r};
    if ~isempty(x)
        plot(x, r*ones(size(x)), 'k.', 'MarkerSize', 8);
    end
end
xlim([t_ms(1) t_ms(end)]);
ylim([0 nrep + 1]);
xlabel('Time from stimulus onset (ms)');
ylabel('Repeat #');
title(sprintf('Location %d', stim_list(stim_select)));
set(gca, 'Box', 'off');

subplot(1,2, 2)
spk1ms = spikeMat{stim_select};                 % nRepeats x 700 (1‑ms bins)
% nrep   = size(spk1ms,1);
% ----- baseline (Hz) from pre + post -----
pre_idx   =  t_ms < 0;
post_idx  =  t_ms >= stim_ms;
base_idx  =  pre_idx | post_idx;       % 200 ms pre + 300 ms post = 500 ms
baseline_rate_rep = sum(spk1ms(:, base_idx), 2) / (pre_ms + post_ms) * 1000; % per‑repeat (Hz)
baseline_rate     = mean(baseline_rate_rep);  % scalar baseline to subtract

% ----- rebin to X ms -----
ncol   = size(spk1ms,2);
keep   = floor(ncol/bin_ms)*bin_ms;    % make divisible by bin-size
spk1ms = spk1ms(:,1:keep);
t_ms1  = t_ms(1:keep);

spkXms = reshape(spk1ms, nrep, bin_ms, []);
spkXms = squeeze(sum(spkXms, 2));      % nRepeats x nBins
t_msX  = t_ms1(1:bin_ms:end);          % left edges of X‑ms bins

mean_rateX    = mean(spkXms, 1) * (1000/bin_ms);   % Hz
mean_rateX_bs = mean_rateX - baseline_rate;        % baseline‑subtracted

% ----- plot -----
yMax = max(mean_rateX_bs);
if yMax <= 0, yMax = eps; end
patch([0 stim_ms stim_ms 0], [0 0 yMax yMax], [0.9 0.9 0.9], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.25); hold on
plot(t_msX, mean_rateX_bs, 'k', 'LineWidth', 1.5); % stairs
yline(0, '--', 'Color', [0.5 0.5 0.5]);

xlabel('Time from stimulus onset (ms)');
ylabel('Firing rate (Hz, baseline-subtracted)');
title(sprintf('Location %d (10‑ms bins, baseline = %.2f Hz)', ...
      stim_list(stim_select), baseline_rate));
box off; axis tight;


