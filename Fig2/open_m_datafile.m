function [ datafile ] = open_m_datafile( monkey_id, file_number, channel_no, tWindow, condition, resolution)

%Open a spike/behavior .m file
%monkey_id is should be 'M71V'
%file_number should be input as a number, like 888
%tWindow is for specifying tmax and tmin, spikes within this window will be summed together
%SPS(97.7kHz) is corrected (only for data_format<5) in the Line 68~70 
%call "TruncateData" "header_reader_mfile" "header_handler" "format_behavior_data"

if file_number < 10
    file_number = ['000' num2str(file_number)];
elseif file_number < 100
    file_number = ['00' num2str(file_number)];
elseif file_number < 1000
    file_number = ['0' num2str(file_number)];
else
    file_number = num2str(file_number);
end

if ~exist('condition','var') || isempty(condition)
    condition = 0;
end
if ~exist('tWindow','var')
    tWindow = [];
end

file_name = [monkey_id file_number];
file_name_m = [file_name '.m'];%to differentiate .m file from .dat file

fid = fopen(file_name_m);
if ~exist(file_name_m,'file')
    datafile = [];
    return
end
if fgetl(fid) == -1
%     datafile = -2;
    fclose(fid);
    datafile = [];
    return;
end

try
    try
        x=eval(file_name);
        data = x.data;
    catch
        x=header_reader_mfile(file_name_m);
        data = x.data;
    end
catch
%     datafile = -2;
    fclose(fid);
    datafile = [];
    return
end

if condition == 1 || condition == 2 %Spike was gained or lost during file
    % file_number
    data = TruncateData(data, channel_no, condition);
    % 'truncating'
end

if ~channel_no
    data(data(:,3) ~= 7,:) = []; %not interested in spikes
end
% spike times need to be corrected b/c circuit does not run at 100kHz but at ~97.7kHz (only for data_format<5)
%thus on each step of the clock the time passed is 1.024 micro seconds instead of 1 micro second
if x.data_format < 5 %non-continuous, not time corrected
    index_nn = find(data(:,4)> -1);
    data(index_nn,4)=data(index_nn,4)*1.024;
else %continuous and already time-corrected. Need to reformat
    %timing information to be trial-based
    if size(data,2) == 6
        behavior_data = data;
    end
    try
        data(:,5:6) = [];
    catch
    end
    try
        data(:,5) = [];
    catch
    end
    %     stimulus_deliveries_indices = find(data(:,3) == 7);
    stimulus_deliveries = data(data(:,3) == 7,:);
    for i = 1:length(stimulus_deliveries)
        %Convert from continuous to trial based timing
        data(data(:,1) == stimulus_deliveries(i,1) & data(:,2) == ...
            stimulus_deliveries(i,2),4) = data(data(:,1) == stimulus_deliveries(i,1) & data(:,2) == ...
            stimulus_deliveries(i,2),4) - (stimulus_deliveries(i,4) - x.pre_stimulus_record_time*1000);
        %I don't think it's necessary to add the delivery indices
        %         for j = 1:6
        %             %add -1 lines (trial markers)
        %             try
        %                 data = [data(1:stimulus_deliveries_indices(i)+(i-1)*6+j-1,:); stimulus_deliveries(i,1)...
        %                     stimulus_deliveries(i,2) j -1; data(stimulus_deliveries_indices(i)+(i-1)*6+j:end,:)];
        %             catch %end of data
        %                 data = [data(1:stimulus_deliveries_indices(i)+(i-1)*6+j-1,:); stimulus_deliveries(i,1)...
        %                     stimulus_deliveries(i,2) j -1];
        %             end
        %         end
    end
end

% Extract Basic Header Parameters
[analysis_type,unit_num,stim_dur,pre_stim,post_stim,hole_num,track_num,unit_depth,] = ...
    header_handler('basic',file_name); %analysis_type,unit_num,stim_dur,pre-stim,post-stim,hole_num,track_num,unit_depth

%there is some serious shit caused by the fact that it's not guaranteed
%that all stimuli in the header will get delivered.
stim = max(max(data(:,1)),length(x.user_stimulus_desc_ch1));
% stimulus_deliveries = data(data(:,3) == 7,:);
% stim = max(stimulus_deliveries(:,1));

if length(stim_dur) > 1 && max(data(:,1)) == length(x.user_stimulus_desc_ch1)%put this in because if two channels, stim dur is returned as twice the length
    stim_dur = stim_dur(1:stim);
elseif max(data(:,1)) == length(x.user_stimulus_desc_ch1) %Weird error in Mod WBE where it doesn't record all of the stim in the header
    stim_dur(1:stim) = stim_dur(1);
end
reps=max(data(:,2));
behav_reps = zeros(1,stim);
if any(data(:,3)==7)
    for i = 1:max(data(:,1))
        behav_reps(i) = sum(data(:,1)==i & data(:,3)==7);
    end
    data(data(:,3)>6,:) = [];
else
    behav_reps = ones(1,stim)*reps;
end

data((data(:,3)~=channel_no | data(:,4)==-1),:)=[]; %get rid of -1 and other channels

if strcmp(resolution, 'ms')
    data(:,4)=round(data(:,4)/1000);     %convert 1000us to 1ms
end

if isempty(tWindow)
    tmin = pre_stim + 15;
    tmax = pre_stim + max(stim_dur) + 20; %original
%     tmax = pre_stim + max(stim_dur) + 40; %used for longer offset firing; CCG @ 2020-05-11
else
    tmin = tWindow{1}(1);
    tmax = tWindow{1}(2);
end
% tmax = pre_stim + max(stim_dur) + post_stim - 50; %users prefer that
% analysis window defaults to stimulus length
if strcmp(resolution, 'ms')
    spont_rate = size(data((data(:,4) < pre_stim & data(:,4) > 0),:),1)/(sum(behav_reps)*pre_stim) * 1000; %spikes/sec
elseif strcmp(resolution, 'us')
    spont_rate = size(data((data(:,4) < pre_stim & data(:,4) > 0),:),1)/(sum(behav_reps)*pre_stim) * 1000000;
end

% Store standard data
datafile.datafile_num = file_number;
datafile.monkey_id = monkey_id;
datafile.date_time = x.datetime;
datafile.analysis_type = analysis_type;
datafile.unit_num = unit_num;
datafile.hole_num = hole_num;
datafile.track_num = track_num;
datafile.unit_depth = unit_depth;
datafile.ch = channel_no;
datafile.stim = stim;
datafile.reps = reps;
datafile.behav_reps = behav_reps;
datafile.stim_dur = stim_dur;
datafile.atten = x.attenuation_ch1;
datafile.stimParams = x.stimulus_ch1;
datafile.pre_stim = pre_stim;
datafile.post_stim = post_stim;
datafile.tmin = tmin;
datafile.tmax = tmax;
datafile.data = data;
datafile.spont_rate = spont_rate;
if ~strcmp(datafile.analysis_type,'RLS')
    [datafile.x_param, datafile.y_param, datafile.dims, datafile.x_name, datafile.y_name]...
        = header_handler('extract_2d',file_name); %gets 2D parameter values
    if iscell(datafile.x_name)
        datafile.x_name = datafile.x_name{1};
    end
    if iscell(datafile.y_name)
        datafile.y_name = datafile.y_name{1};
    end
else
    datafile.x_param = 1:24;
    datafile.y_param = 0;
    datafile.dims = 1;
    datafile.x_name = 'RLS_stim';
    datafile.y_name = 'stim_#';
    indexLoc = 0;
    %need to automatically figure out which files are which:
    %Unique vs. full: number of stim (51 vs. 255)
    %Seeds: single: all same seed value (look for a large number of
    % the same value in the header stimulus information), fixed
    % seed: either "index" parameter doesn't exist, or it does
    % exist and every channel starts on the same index in every
    % timulus, unique: index parameter exists and each channel
    % starts at a different index in each seed.
    for i = 1:numel(x.stimulus_tags_ch1)
        if strcmp(' Indexes',x.stimulus_tags_ch1{i})
            indexLoc = i;
            %found indexes, go throuh indexes of first stim to see whether
            %fixed (including single) or unique seeds.
            if x.stimulus_ch1(1,indexLoc) == x.stimulus_ch1(2,indexLoc)
                datafile.rlsSeed = 'single';
            elseif x.stimulus_ch1(1,indexLoc) == x.stimulus_ch1(1,indexLoc+1)
                datafile.rlsSeed = 'fixed';
            else
                datafile.rlsSeed = 'unique';
            end
            break
        end
    end
    if ~indexLoc
        datafile.rlsSeed = 'fixed';
    end
    if size(x.stimulus_ch1,1) > 51
        datafile.rlsStim = 'unique';
    else
        datafile.rlsStim = 'not unique';
    end
end
if exist('behavior_data','var')
    datafile.behavior_data = behavior_data;
    datafile.behavior_tags = x.behavior_tags;
    datafile.behavior_params = x.behavior.params;
    datafile.behavior_type = x.behavior_type;
    datafile.sham_percentage = x.sham_percentage;
    datafile.behav_iti_min = x.behav_iti_min;
    datafile.behav_iti_max = x.behav_iti_max;
    datafile.blackout_dur = x.blackout_dur;
    datafile.air_puff = x.air_puff;
    if channel_no %looking for spike data
%         try
            [datafile.behavior_data_formatted, datafile.behavior_reps_formatted] = format_behavior_data(datafile);
%         catch err
%             if strcmp(err.identifier,'MATLAB:badsubscript')
%                 disp('behavior datafile subscript error')
%                 datafile = [];
%                 return
%             else
%                 rethrow(err)
%             end
%         end
    end
end

fclose(fid);

end


