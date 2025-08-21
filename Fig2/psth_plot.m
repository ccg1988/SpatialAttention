function [t_msX, mean_rateX_bs, baseline_rate] = psth_plot(spikeMat, stim_select, bin_ms, pre_ms, stim_ms, post_ms, t_ms)

spk1ms = spikeMat{stim_select};                 % nRepeats x 700 (1‑ms bins)
nrep   = size(spk1ms,1);
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
