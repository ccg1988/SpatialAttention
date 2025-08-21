function [spikeTimes, spikeMat, pre_ms, stim_ms, post_ms, t_ms, stim_list] = spike_extract(file_name, CH)

stim_col  = 1;   % stimulus ID
label_col = 3;   % event code: 4=spike, 7=stimulus onset, 8/9=lick, 10=reward
time_col  = 4;   % time in microseconds

% file_name='M9X0831';
data_raw=str2num(file_name);
data=data_raw.data;
S=data_raw.stimulus_ch1; 
time_ms = double(data(:, time_col)) / 1000; % Convert us to ms
is_spike = data(:, label_col) == CH ; % Event masks
is_stim  = data(:, label_col) == 7 ;
pre_ms    = data_raw.pre_stimulus_record_time;           % prestimulus window
stim_ms   = S(1,5);           % stimulus duration
post_ms   = data_raw.post_stimulus_record_time;           % poststimulus window
t_ms   = -pre_ms : (stim_ms + post_ms - 1);   % e.g., -200 ... 499  (700 bins); 1 ms resolution
edges  = [t_ms, t_ms(end) + 1];               % histcounts edges
stim_list = unique(S(:,2), 'stable') ;   % five stimuli and first is always "background"
spikeTimes = cell(numel(stim_list), 1);  % cell{stim}{repeat} = vector of spike times (ms, relative to onset)
spikeMat   = cell(numel(stim_list), 1);  % cell{stim} = [nRepeats x 700] counts per ms
for si = 1:numel(stim_list)
    sid = stim_list(si);
    ids = S(S(:,2) == sid, 1);                      % map type (8 vs 3-7-15-17) -> its 11 true IDs
    onset_idx = is_stim & ismember(data(:, stim_col), ids);
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