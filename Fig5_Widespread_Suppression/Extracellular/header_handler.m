function varargout = header_handler(control_str,file_name,varargin)

%This function extracts relevant information from the data file headers (.m & .ad2 file).
%.ad2 file: Line 8   to 86
%.m file:    Line 88 to end  "quick, basic, extract_2d (very complex), extract_ici, stim_description, rss, rls, data_format"
%Line 107: when data.format<5 (SPS=97.7), do noting. However, it is corrected in the "open_m_datafile" 


%If looking for .ad2 header, the whole header will be passed in instead of the file name
if strcmp(lower(control_str),'ad2')
    header = file_name;
    %Data file format version
    %AD2 channel on flags
    %AD2 channel labels
    i = 2;
    j = 0;

    % Search forwards
    while  j < 1 %looking for 1 line
        s = header{i};
        if j == 0
            found = findstr(s,'Data file format');
            if found
                varargout{1} = get_numbers(s);
            end
        end
        if found %look for the next thing
            j = j+1;
        end

        i = i+1;
        if i>length(header)
            i=2;
            varargout{1} = 1;
            break
        end
    end %while

    i = length(header)-1;
    j = 0;

    % Search backwards
    while  j < 2 %looking for 2 lines
        s = header{i};
        if j == 0
            found = findstr(s, 'AD2 Channel Labels');
            if found
                [t,r] = strtok(s,'=');
                [t,r] = strtok(r);
                ad2_labels = r;
            end
        elseif j == 1
            found = findstr(s, 'AD2 Channel On Flags');
            if found
                [t,r] = strtok(s,'=');
                [t,r] = strtok(r);
                on_ch = str2num(r);
                ad2_on_ch = find(on_ch);

                clear t
                for k = 1:length(on_ch)
                    [t{k},ad2_labels] = strtok(ad2_labels);
                end

                ad2_labels = [];
                for k = 1:length(ad2_on_ch)
                    ad2_labels = strvcat(ad2_labels,t{ad2_on_ch(k)});
                end

                varargout{2} = ad2_labels;
                varargout{3} = ad2_on_ch;
            end
        else
            %nothing left to look for
        end
        if found %look for the next thing
            j = j+1;
        end
        i = i-1;
        if i<2
            i=length(header)-1;
            varargout{2} = 'Spikes';
            varargout{3} = 1;
            break
        end
    end %while
    return
end
%************************************************for .m file******************************************************
try  %% only use header_reader_mfile if necessary - i.e. updating zero while running stimuli.
    x = eval(file_name);
catch
    x=header_reader_mfile([file_name '.m']);
end

switch lower(control_str)
    case 'quick'
        varargout{1} = x.analysis_type;
        varargout{2} = x.unit_number;
        varargout{3} = x.nstim;
        varargout{4} = mean(x.stimulus_ch1(:,5));
        varargout{5} = max(x.data(:,2));
        
        data = x.data;
        
        if x.data_format < 5 %non-continuous, not time corrected
            %do nothing
        else %continuous and already time-corrected. Need to reformat
            %timing information to be trial-based
            stimulus_deliveries_indices = find(data(:,3) == 7);
            stimulus_deliveries = data(data(:,3) == 7,:);
            for i = 1:length(stimulus_deliveries)
                data(data(:,1) == stimulus_deliveries(i,1) & data(:,2) == ...
                    stimulus_deliveries(i,2),4) = data(data(:,1) == stimulus_deliveries(i,1) & data(:,2) == ...
                    stimulus_deliveries(i,2),4) - (stimulus_deliveries(i,4) - x.pre_stimulus_record_time*1000);
                for j = 1:6
                    try
                        data = [data(1:stimulus_deliveries_indices(i)+(i-1)*6+j-1,:); stimulus_deliveries(i,1)...
                            stimulus_deliveries(i,2) j -1; data(stimulus_deliveries_indices(i)+(i-1)*6+j:end,:)];
                    catch
                        data = [data(1:stimulus_deliveries_indices(i)+(i-1)*6+j-1,:); stimulus_deliveries(i,1)...
                            stimulus_deliveries(i,2) j -1];
                    end
                end
            end
            data(data(:,3)>6,:) = [];
        end
        
        stim1_chs = data(find(data(:,1)==1 & data(:,2)==1) & data(:,3) < 6,3);
        chs(1) = stim1_chs(1);
        ch_dist = stim1_chs - chs(1);
        ch_dist = ch_dist(find(ch_dist(2:end) - ch_dist(1:end-1) > 0)+1);
        chs = chs(1) + [0 ch_dist'];
        
        varargout{6} = chs;
        
    case 'basic'
        varargout{1} = x.analysis_type;
        varargout{2} = x.unit_number;
        varargout{3} = x.stimulus_ch1(:,5);
        varargout{4} = x.pre_stimulus_record_time;
        varargout{5} = x.post_stimulus_record_time;
        varargout{6} = x.hole_number;
        varargout{7} = x.track_number;
        varargout{8} = x.unit_depth - x.starting_depth;
        
    case 'extract_2d' %try to find 2 independent variables in the stimulus portion of the header
        analysis_code = x.analysis_code;
        if analysis_code >= 700 & analysis_code < 800
            is_rss = 1; %RSS stim headers use multiple lines per stim
        else
            is_rss = 0;
        end
        
        stim_params = x.stimulus_ch1;
        
        % Get rid of [stim | speaker | reps | duration] since these parameters vary all of the time
        stim_params(:,[1,4:5])=[];
        
        %	d_stim_params = stim_params(1,:) - stim_params(end,:); %look at difference of first and last stim
        %	d_index_outer = find(d_stim_params~=0); % d_index is the indices of the changing parameters
        
        %% Determine Dimensions - Begin
        full_d_stim_params = stim_params(2:end,:) - stim_params(1:(end-1),:); %subtract matrix from itself
        full_d_stim_params(end+1,:) = stim_params(1,:) - stim_params(end,:); %wrap around
        
        full_d_stim_params(find(full_d_stim_params)) = 1; %1's indicate a change of param from previous trial
        
        num_changes = sum(full_d_stim_params,1); %a row vector
        
        % Get rid of two or more simultaneously changing columns
        j = find(num_changes);
        for i = 1:length(j)
            temp_num_changes = num_changes - num_changes(j(i)); %subtract the current # of changes
            temp_num_changes(j(i)) = num_changes(j(i)); %keep at least the 1st element with that # of changes
            num_changes(find(temp_num_changes==0)) = 0; %set extra columns with same # of changes to 0
        end
        
        % If one parameter accounts for all changes, but other parameters
        % vary in a non-multiplicative way
        for i = 1:length(num_changes)
            if num_changes(i) ~= 0
                if size(stim_params,1)/num_changes(i) ~= round(size(stim_params,1)/num_changes(i)) %if not an integer multiple
                    num_changes(i) = 0;
                end
            end
        end
        
        dims = length(find(num_changes));
        total = size(stim_params,1);
        
        num_changes((num_changes==0)) = total + 1; %insures that counts of 0 are not the min
        
        % look for outermost changing parameters
        if dims >= 2 %pick two params
            [v,i] = min(num_changes); %outermost param
            x_ind = i;
            jump = total/num_changes(i);
            x = stim_params(jump:jump:end,x_ind);
            
            num_changes(i) = total + 1; %insures that this index will not be the min again
            
            [v,i] = min(num_changes); %2nd outermost param
            y_ind = i;
            jump = total/num_changes(i); %both parameters overcount by the length of the outermost parameter
            total = total/length(x); %just look at a block where outermost parameter does not change, since that will yield one full cycle of 2nd outer param
            y = stim_params(jump:jump:total,y_ind);
        elseif dims == 1 %only one varying param
            [v,i] = min(num_changes);
            x_ind = i;
            x = stim_params(:,x_ind);
            y = 0;
            y_ind = 0;
        elseif dims == 0 %just use stim # as param
            x = [1:size(stim_params,1)];
            y = 0;
            x_ind = 0;
            y_ind = 0;
        end
        %% Determine Dimensions - End
        
        x_label = get_independent_variable_label(x_ind,analysis_code);
        y_label = get_independent_variable_label(y_ind,analysis_code);
        
        varargout{1} = x;
        varargout{2} = y;
        varargout{3} = dims;
        varargout{4} = x_label;
        varargout{5} = y_label;
        if x_ind ~= 0
            varargout{6} = stim_params(:,x_ind);
        else
            varargout{6} = [];
        end
        
        if y_ind ~= 0
            varargout{7} = stim_params(:,y_ind);
        else
            varargout{7} = [];
        end
    case 'extract_ici'
        
        % figure out whether there is actually a temporally modulated
        
        % look for an ici or an fmod and find out column #
        temp_mod_str = {'ici','fmod','modulation'};
        temp_mod_units = {0,1,1}; % 0 means in ms 1 means in Hz
        mod_flag = 0;
        for i1 = 1:size(temp_mod_str,2)
            i_mod = strmatch(temp_mod_str{i1},strtrim(lower(x.stimulus_tags_ch1)));
            if ~isempty(i_mod)
                mod_flag = 1;
                col_mod = i_mod;
                mod_units = temp_mod_units{i1};
            end
        end
        if mod_flag == 0 % stop if no modulation has been found
            varargout{1} = [];
            return
        end
        
        % determine ICIs
        i_mod = [0:1:x.nstim-1]*size(x.stimulus_tags_ch1,2)+col_mod
        ici = x.stimulus_ch1(i_mod);
        if mod_units == 1
            ici = 1./ici.*1000; %convert from Hz to ms
        end
        varargout{1} = ici;
    case 'stim_description'
        stim_num = round(varargin{1})+1;
        
        stim_desc{stim_num,1} = x.user_stimulus_desc_ch1{stim_num};
        
        if ~isempty(x.stimulus_ch2)
            stim_desc{stim_num,2} = x.user_stimulus_desc_ch2{stim_num};
        else
            stim_desc{stim_num,2} = [];
        end
        varargout{1} = stim_desc;
    case 'rss'
        level_matrix = [];
        for i = 1:x.nstim;
            stim_str = x.user_stimulus_desc_ch1{i};
            level_matrix(i,:) = str2num(stim_str(findstr(stim_str,'StimLevs=')+10:end));
        end
        varargout{1} = level_matrix;
        
        SD_index = findstr(stim_str,'LevSD(dB)=')+10;
        SD = str2num(stim_str(SD_index:SD_index+1));
        varargout{2} = SD;
        
        c_freq_index = findstr(stim_str,'Freq(Hz)=')+9;
        c_freq = str2num(stim_str(c_freq_index:c_freq_index+3))/1000;
        BW_index = findstr(stim_str,'BW(oct)=')+8;
        BW = str2num(stim_str(BW_index));
        
        freq = 2.^([log2(c_freq)-BW/2:BW/(size(level_matrix,2)-1):log2(c_freq)+BW/2]);
        varargout{3} = freq;
        
    case 'rls'
        level_matrix = [];
        for i = 1:x.nstim;
            stim_str = x.user_stimulus_desc_ch1{i};
            level_matrix(i,:) = str2num(stim_str(findstr(stim_str,'Stimulus Levels')+16:findstr(stim_str,'Seeds')-1));
        end
        varargout{1} = level_matrix;
        
        SD_index = findstr(stim_str,'Standard Deviation')+19;
        SD = str2num(stim_str(SD_index:SD_index+1));
        varargout{2} = SD;
        
    case 'data_format'
        varagout{1} = x.data_format;
end %switch

end

