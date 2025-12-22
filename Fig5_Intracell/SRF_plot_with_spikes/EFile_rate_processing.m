function rate_structure = EFile_rate_processing(x,ch_number,onset_time,offset_time)

% input x: original data structure x
% input onset/offset time: time window to calculate firing rate
% output: rate_structure


data = x.data;  

pre_stim = x.pre_stimulus_record_time;  % in mS
post_stim = x.post_stimulus_record_time; % in mS
stim_dur = x.stimulus_ch1(1,5); % in mS
S_raw=x.stimulus_ch1; S=S_raw(:,6); atten=S_raw(:,3);
Rep = S_raw(:,4); Rep = max(Rep);

tmin = pre_stim + onset_time;
tmax = pre_stim + stim_dur + offset_time;
    
stim = max(data(:,1));
for i = 1:stim
    reps(i) = max(data(data(:,1) == i,2));
end

spkr_number_vector = [];

if strcmp(x.trial_complete,'true') ==  0
    % delete unfinished repetition's data
    data(data(:,2)>min(reps),:)=[];
end
% data(data(:,3)~=ch_number |data(:,4)==-1,:)=[]; %get rid of -1 and other channels
data(data(:,3)~=ch_number | data(:,4)<=0,:)=[]; % for both extra(-1) and intra(-40) recording in Chamber 1
data(:,4)=round(data(:,4)/1000); %convert to ms
for i = 1:stim
    spkr_number_vector(i) = x.stimulus_ch1(i,2);
    for j = 1:min(reps)
        spont_train = data(data(:,1)==i & data(:,2)==j & data(:,4) <= pre_stim, 4);   %in mS 
%         spont_train = data(data(:,1)==i & data(:,2)==j & data(:,4) <= pre_stim & data(:,4) > 0, 4);   %CCG 20200710
        spike_train = data(data(:,1)==i & data(:,2)==j & data(:,4)>=tmin & data(:,4) <= tmax,4);   %in mS 
        spont_rates(i,j) = length(spont_train)/pre_stim; % spikes/ms
               
        driven_spcounts(i,j) = length(spike_train); % spike counts
        rel_driven_spcounts(i,j) = length(spike_train) - length(spont_train)/pre_stim*(tmax-tmin);
        clear spont_train spike_train;
    end
end

rate_structure.spkr_number_vector = spkr_number_vector;
rate_structure.spont_rates = spont_rates;  %spikes/ms
rate_structure.driven_spcounts = driven_spcounts;
rate_structure.rel_driven_spcounts = rel_driven_spcounts;
rate_structure.analysis_onset_time = onset_time; %ms
rate_structure.analysis_offset_time = offset_time; %ms
rate_structure.analysis_window_length = tmax-tmin;%ms
rate_structure.stim_dur = stim_dur;
rate_structure.spont_analysis_window = pre_stim; %ms
rate_structure.type = x.analysis_type; % 'WB-noise' or 'bandwidth'
rate_structure.type = S; % the stimuli used (frequency or location)
rate_structure.atten = atten; 
rate_structure.rep = Rep; 
rate_structure.side = x.hemisphere;