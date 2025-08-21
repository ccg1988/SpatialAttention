function [ datafile, best_location, number_of_peaks, tuning_area, tuning_vector_magnitude,...
    rate_interpolation, srf_hfig] = analyze_srf_beta( varargin )

% Simplified by CCG @ 2020-12-03
% This function does NOT call other functions, but need loading speaker information

% Functionality for inputing individual input requirements is deprecated.

% The functions of divide_by_zero have two working modes: 
% In 1st mode, a datafile structure is input, which will should have all the necessary information. 
% In 2nd mode, the needed information is input individually.
% For individual input args, the function call should look like this:
% [best_loc num_peaks tuning_width tuning_vector rate_interp ] = analyze_SRF_beta(rates,0,spont_rate,unit_number,plot_on);

%If the rates are given directly to the function, no spont rate subtraction is made, you can do this yourself. Otherwise, rates are spont-rate corrected

%tuning vector is [left(negative)/right(positive) front(positive)/rear(negative) up(positive)/down(negative)]

%% Get data

need_data = 1;

if isstruct(varargin{1}) %passed in datfile struct
    datafile = varargin{1};
    params_used = 1;
    extra_params = length(varargin) - params_used;
elseif ischar(varargin{1}) %passed in monkey ID and file number/channel number
    datafile = open_m_datafile(varargin{1},varargin{2}(1), varargin{2}(2));
    params_used = 2;
    extra_params = length(varargin) - params_used;
else %passed in data directly
    rates = varargin{1} - varargin{2}; %subtract spontaneous rate
    if size(rates,2) == 1
        rates = rates';
    end
    spont_rate = varargin{2};
    stdev_of_mean = zeros(1,length(rates)); %can't put in error bars
    params_used = 2;
    extra_params = length(varargin) - params_used;
    need_data = 0;
end

if need_data
    unit_number = datafile.unit_num;        %original is number
    file_number = datafile.datafile_num;   %original is number
    ch = datafile.ch;
    if isfield(datafile,'rls_weights')
        y_axis_label = 'Spikes/Second/dB';
        figure_title = ['RLS Weights, Unit ' num2str(unit_number) ', Channel ' num2str(ch) ', File Number  ' num2str(file_number)];
        rates = datafile.rls_weights;
        spont_rate = -1*min(rates);
        stdev_of_mean = 0;
    else
        y_axis_label = 'Spikes/Second';
        figure_title = ['File: ',num2str(file_number),'  U: ' num2str(unit_number) '  Ch: ' num2str(ch),...
            '  R: ',num2str(size(datafile.rep_rates,1)),'  Att: ',num2str(datafile.atten)];
        if isfield(datafile,'rates')
            rates = datafile.rates - datafile.spont_rate;
            stdev_of_mean = datafile.stdev_of_mean;
        else
            [datafile] = analyze_rates(datafile);
            rates = datafile.rates - datafile.spont_rate;
            stdev_of_mean = datafile.stdev_of_mean;
            y_axis_label = 'Spikes/Second';
        end
        spont_rate = datafile.spont_rate;
    end
end

if extra_params && ~isempty(varargin{params_used+1})
    plot_on = varargin{params_used+1};
else
    plot_on = 0;
end
if extra_params > 1 && ~isempty(varargin{params_used+2})
    orientation = varargin{params_used+2};
elseif exist('datafile','var') && isfield(datafile,'orientation')
    switch datafile.orientation
        case 'facing front'
            orientation = 1;
        case 'facing left'
            orientation = 2;
        case 'facing right'
            orientation = 3;
        case 'facing back'
            orientation = 4;
    end
else
    orientation = 1;
end
human_name='Marmoset';
if extra_params > 2 && ~isempty(varargin{params_used+3})
    human_name=varargin{params_used+3};%CCG
    behav_speakers = [];%CCG
%     behav_speakers = varargin{params_used+3};%original
else
    behav_speakers = [];
end
if extra_params > 3 && ~isempty(varargin{params_used+4})
    y_axis_label = varargin{params_used + 4};
elseif ~exist('y_axis_label','var')
    y_axis_label = ' ';
end
if extra_params > 4 && ~isempty(varargin{params_used+5})
    figure_title = varargin{params_used + 5};
elseif ~exist('figure_title','var')
    figure_title = ' ';
end
if extra_params > 5 && ~isempty(varargin{params_used+6})
    threshold_type = varargin{params_used + 6};
else
    threshold_type = 1;
end


%% Initialize outputs
best_location = [];
number_of_peaks = [];
tuning_area = [];
tuning_vector = [];
rate_interpolation = [];
srf_hfig = [];

%% Setup spatial grid
if orientation == 1
    srf = load('SRF.mat');
else
    srf = load(['SRF' num2str(orientation) '.mat']);%SRF2.mat or SRF3.mat
end
if strcmp (human_name,'YYW')
   srf=load('SRF_YYW.mat');
end
if strcmp (human_name,'ER')
   srf=load('speakers1_32_ER.mat');
end
speaker_locations = srf.speakers;
speaker_number=size(speaker_locations,1);
spatial_resolution = 5;
longitude_limit = 180-spatial_resolution/2;
latitude_limit = 90-spatial_resolution/2;
nlongitude_cells = 360/spatial_resolution;
nlatitude_cells = 180/spatial_resolution;

%% Do calculations based on number of locations
if strcmp (human_name,'T')
    rates = rates/max(rates);
end


if length(rates) == speaker_number %normal 24 or 32 locations SRF
    
    longitude = repmat(-1*longitude_limit:spatial_resolution:longitude_limit,nlatitude_cells,1);
    latitude = repmat(-1*latitude_limit:spatial_resolution:latitude_limit,nlongitude_cells,1)';
    
    latitude_up = latitude(:,1) + spatial_resolution/2;
    latitude_down = latitude(:,1) - spatial_resolution/2;
    longitude_left = longitude(:,1) + spatial_resolution/2;
    longitude_right = longitude(:,1) - spatial_resolution/2;

    interp_cell_areas = areaquad(latitude_down,longitude_right,latitude_up,longitude_left);
    interp_cell_areas = repmat(interp_cell_areas,1,nlongitude_cells);
    
    distances = zeros(speaker_number, nlatitude_cells, nlongitude_cells);
    for i = 1:speaker_number %faster to do it this way instead of element by element. Takes .01 - .3 seconds
        distances(i,:,:) = reshape(distance(speaker_locations(i,[2 1]), [latitude(:) longitude(:)]),nlatitude_cells,nlongitude_cells);
    end
      
    rate_interpolation = zeros(nlatitude_cells,nlongitude_cells);
    for i = 1:nlatitude_cells
        for j = 1:nlongitude_cells
            [dist, index] = sort(distances(:,i,j));
            weights = (1./dist(1:2)).^2./(sum(1./dist(1:2).^2));
            rate_interpolation(i,j) = sum(weights'.*rates(index(1:2)));
        end
    end
    
    if strcmp (human_name,'T')
%     rate_interpolation = flipud(rate_interpolation);
    rate_interpolation=rate_interpolation/max(rate_interpolation,[],'all');
    C_x=36; C_y=18;
    [~,peak_x] = max( (max(rate_interpolation,[],1)) );
    [~,peak_y] = max( (max(rate_interpolation,[],2)) );
    T_x=C_x-peak_x;
    T_y=C_y-peak_y;
    temp_x = circshift (rate_interpolation ,T_x, 2 ); %2 for X-demension
    rate_interpolation=imtranslate(temp_x,[0, T_y]);
    end
    
    [rate_max, best_location] = max(rates);
    %The new threshold is half the maximum firing rate (raw, not spont corrected). 
    %The fussing with spont here is because the rates in this function are spont subtracted.
    if threshold_type == 1
        threshold = .5*(rate_max+spont_rate)-spont_rate;
    elseif threshold_type == 2
        threshold = .5*(rate_max);
    end
    
    threshold_contour = contourc(rate_interpolation,[threshold threshold]);
    index = 1;
    edge_touches = [0 0];

    %%
    %New Way (Marcus)
    number_of_contours = 0;
    final_contour = 0;
    to_clear_index = [];
    i = 1;
    while final_contour == 0 && ~isempty(threshold_contour)
        contour_length = threshold_contour(2,index);
        if contour_length > 10
            threshold_contour(:,index) = NaN;
            %% test if contour encloses an area above or below threshold
            raw_contours{i} = threshold_contour(:,index+1:index+contour_length);
            if any(raw_contours{i}(1,:) == 1)
                [~,contour_index] = max(raw_contours{i}(1,:));
                rate_interp_index(2) = floor(raw_contours{i}(1,contour_index));
                rate_interp_index(1) = round(raw_contours{i}(2,contour_index));
            else
                [~, contour_index] = min(raw_contours{i}(1,:));
                rate_interp_index(2) = ceil(raw_contours{i}(1,contour_index));
                rate_interp_index(1) = round(raw_contours{i}(2,contour_index));
            end
            threshold_test = rate_interpolation(rate_interp_index(1), rate_interp_index(2));
            if threshold_test > threshold
                number_of_contours = number_of_contours + 1;
                %just x values for scaled_contours
                scaled_contours{i} = threshold_contour(1,index+1:index+contour_length)* spatial_resolution-(180+spatial_resolution/2);
                scaled_contours{i}(scaled_contours{i}(1,:) == max(max(longitude))) = 180;
                scaled_contours{i}(scaled_contours{i}(1,:) == min(min(longitude))) = -180;
                scaled_contours{i} = scaled_contours{i}(1,:)/180;
                %look for contours that go from one side to the other so we don't
                %overestimate the number of peaks in the SRF
                %way to check:
                %1: Both edges touching same side
                %2: only one edge is touching
                if (sum(abs(scaled_contours{i}([1 end])) ==  1) == 2 && scaled_contours{i}(1) == scaled_contours{i}(end)) || ...
                        sum(abs(scaled_contours{i}([1 end])) ==  1) == 1 %%touched, wraps around
                    if any(scaled_contours{i}([1 end]) == -1) %-180
                        edge_touches(1) = edge_touches(1) + 1;
                    else % 180
                        edge_touches(2) = edge_touches(2) + 1;
                    end
                end
            end
        else
            to_clear_index = [to_clear_index [index:index + contour_length]];
        end
        index = index + contour_length + 1;
        if index >= size(threshold_contour,2)
            final_contour = 1;
        end
        i = i + 1;
    end
    threshold_contour(:,to_clear_index) = [];
    
    % threshold_contourScaled = threshold_contour./[nlongitude_cells nlatitude_cells];
    threshold_contour(1,:) = threshold_contour(1,:)*spatial_resolution-(180+spatial_resolution/2);
    threshold_contour(2,:) = threshold_contour(2,:)*spatial_resolution-(90+spatial_resolution/2);
    threshold_contour(1,threshold_contour(1,:) == max(max(longitude))) = 180;
    threshold_contour(1,threshold_contour(1,:) == min(min(longitude))) = -180;
    threshold_contour(2,threshold_contour(2,:) == max(max(latitude))) = 90;
    threshold_contour(2,threshold_contour(2,:) == min(min(latitude))) = -90;
    
    % edge_touches = sum(any(threshold_contourNorm == 1,1));
    
    number_of_peaks = number_of_contours - max(edge_touches);
    if ~isempty(threshold_contour)
        number_of_peaks = max(number_of_peaks,1);
    end
    
%     if number_of_peaks > 4
%         plot_on = 2;
%     end
    
    %just use vector (1-D) indexing for the matrix
%     cells_above_threshold = find(rate_interpolation > threshold);
    
%     latitude_up = latitude(cells_above_threshold) + spatial_resolution/2;
%     latitude_down = latitude(cells_above_threshold) - spatial_resolution/2;
%     longitude_left = longitude(cells_above_threshold) + spatial_resolution/2;
%     longitude_right = longitude(cells_above_threshold) - spatial_resolution/2;
    tuning_area = sum(interp_cell_areas(rate_interpolation > threshold));
    
    
    %*********************************Tuning vector and angles*********************************
    %tuning_vector = mean([rates' rates' rates'].*speaker_locations(:,3:5))/mean(rates);---Wrong
    % old tuning vector is biased by speaker distribution (analyze_SRF, L: 231)
    % tuning_vector = (rates*speaker_locations(:,3:5))/sum(rates-min(rates));
    
    interp_vectors = zeros(nlatitude_cells*nlongitude_cells,3);
    [interp_vectors(:,2), interp_vectors(:,1), interp_vectors(:,3)] = ...
        sph2cart(pi*longitude(:)/180,pi*latitude(:)/180, ...
        (rate_interpolation(:)+spont_rate).*interp_cell_areas(:));
    tuning_vector = 1.066*sum(interp_vectors)/sum((rate_interpolation(:)+spont_rate).*interp_cell_areas(:));
    %1.066 is an adjustment so that a response to only one location has a magnitude of 1;
    [~,~,tuning_vector_magnitude] = cart2sph(tuning_vector(2),tuning_vector(1),tuning_vector(3));
    tuning_vector_magnitude = min(tuning_vector_magnitude,1);
    %%
    [lateral_centroid, best_lateral_angle, yz_centroid]= horizontal_pole_transform(rates, 0);% 0="Off", 1="On"
    %%
%     if tuning_area < .3 && tuning_vector_magnitude < .3
%         plot_on = 1;
%     end
%**************************************************************Plot SRF*****==1/2/3******************************************
    if plot_on  
        if plot_on == 1
            srf_hfig = figure('name',figure_title,'position',[30 250 600 325]);
        elseif plot_on == 2
            srf_hfig = figure('name',figure_title,'position',[30 250 525 525]);
        elseif plot_on == 3
            srf_hfig = figure('name',figure_title,'position',[30 250 820 525]);
        end
        
        spatialR = georasterref('RasterSize', size(rate_interpolation), 'Latlim', [-90 90], 'Lonlim', [-180 180]); 
        %spatialref.GeoRasterReference object
        ax = axesm('MapProjection', 'mollweid');%fournier is used by Evan's paper
%         setm(ax,'Grid','on','Frame','off','MeridianLabel','On','MlabelParallel',10, ...
%             'ParallelLabel','On','LabelFormat','none');
        setm(ax,'Grid','on','Frame','off','MlabelParallel',10,'LabelFormat','none');
        axis off
        
        if plot_on == 3
            set(gca,'position',[.02 .5 .64 .5]);
        elseif plot_on == 2
            set(gca,'position',[.02 .5 1 .5]);
        else
            set(gca,'position',[0.02 0 1 1]);
        end
        %Show the inter-polated spike rates
        geoshow(rate_interpolation, spatialR,'DisplayType', 'texturemap'); 
        colormap(jet)
        absolute_maximum = max(abs(rates));
        if length(unique(rates)) > 1 %Avoiding error
            set(gca,'Clim',[-1*absolute_maximum absolute_maximum]);
        end
        
        hemisphere_line = [-90*ones(181,1) (-90:90)'; 90*ones(181,1) (90:-1:-90)'];
        %%%Show the white-circle hemisphere line
        geoshow(hemisphere_line(:,2), hemisphere_line(:,1),'color',[1 1 1],'linewidth',2);
        
        location_numbers = 1:speaker_number;
        locations_index = zeros(1,speaker_number);
        if ~isempty(behav_speakers)
            background_location = locations_index;
            background_location(behav_speakers(1)) = 1;
            number_of_backgrounds = sum(background_location);
            target_locations = locations_index;
            target_locations(behav_speakers(2:end)) = 1;
            number_of_targets = sum(target_locations);
            passive_locations = ~(background_location | target_locations);
            number_of_passive = sum(passive_locations);
        else
            passive_locations = ones(1,speaker_number);
            number_of_passive = speaker_number;
        end
        
        locations = struct('Geometry','Point', ...
            'lat',mat2cell(speaker_locations(~~passive_locations,2),ones(1,number_of_passive)), ...
            'long',mat2cell(speaker_locations(~~passive_locations,1), ones(1,number_of_passive)), ...
            'Name',cellstr(num2str(location_numbers(~~passive_locations)')));
        %Show speaker locations
        geoshow(locations,'Marker','o', 'MarkerFaceColor',[0 0 0], 'MarkerEdgeColor','k', 'MarkerSize',7);
        
        if  ~isempty(behav_speakers)
            background_geostruct = struct('Geometry','Point', ...
                'lat', mat2cell(speaker_locations(~~background_location,2),ones(1,number_of_backgrounds)),...
                'long',mat2cell(speaker_locations(~~background_location,1), ones(1,number_of_backgrounds)), ...
                'Name',cellstr(num2str(location_numbers(~~background_location)')));
            geoshow(background_geostruct,'Marker','o', 'MarkerFaceColor', [.5 .5 .5],'MarkerEdgeColor','k','MarkerSize',20);
            %         textm([background_geostruct(:).Lat]*.8,[background_geostruct(:).Lon],'B','FontSize',12)
            
            targets_geostruct = struct('Geometry','Point', ...
                'lat', mat2cell(speaker_locations(~~target_locations,2),ones(1,number_of_targets)),...
                'long',mat2cell(speaker_locations(~~target_locations,1),ones(1,number_of_targets)), ...
                'Name',cellstr(num2str(location_numbers(~~target_locations)')));
            geoshow(targets_geostruct,'Marker','o', 'MarkerFaceColor',[1 1 1], 'MarkerEdgeColor','k', 'MarkerSize',20);
            %         textm([targets_geostruct(:).Lat]*.8,[targets_geostruct(:).Lon],'T','FontSize',12)
        end
        %%
        if plot_on == 1
            text(.4, -1.05, {['Lateral Centroid: ' num2str(lateral_centroid,2)], ['Vertical Centroid: ' num2str(yz_centroid,3)]});
            text(.4, 1.05, {['Tuning Area: ' num2str(tuning_area,2)], ['Vector Magnitude: ' num2str(tuning_vector_magnitude,2)]});
        end    
        if length(unique(rates)) > 1
            %Show threshold contour/boundary
            geoshow(threshold_contour(2,:),threshold_contour(1,:),'color','k','linewidth',2);
            %Show Tuning Vector (white-color spot) (size dependent on |TV|)
            [THETA,PHI,~] = cart2sph(tuning_vector(2),tuning_vector(1),tuning_vector(3));
            tuning_vector_geostruct = struct('Geometry','Point', 'lat',rad2deg(PHI), 'long', rad2deg(THETA), 'Name','TV');
            geoshow(tuning_vector_geostruct,'Marker','o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], ... 
                'MarkerSize',20*sqrt(tuning_vector_magnitude));
            textm(tuning_vector_geostruct.lat, tuning_vector_geostruct.long-9*tuning_vector_magnitude, ...
                tuning_vector_geostruct.Name, 'FontSize',14*tuning_vector_magnitude);

            %Show scale-bar of spike/second
           colorbar('Ylim',[min(rates) max(rates)],'location','westoutside');
           if strcmp (human_name,'T')
               colorbar('Ylim',[-1 1],'location','westoutside');
           end
            x_position = min(get(gca,'xlim'))-0.5;%shift "y_axis_label" to the left
            y_position = min(get(gca,'ylim'));
            text(x_position,y_position,y_axis_label,  'HorizontalAlignment','Left', 'VerticalAlignment','Top');       
        end       
        %%
        if plot_on > 1
            text(.7, -1.6, {['Tuning Area: ' num2str(tuning_area,2)], ['Vector Magnitude: ' num2str(tuning_vector_magnitude,2)]});
            if plot_on == 2
                hrates = subplot(2,1,2);
                hrates_position = get(hrates,'Position');
            elseif plot_on == 3
                srf_hfig(3) = subplot(1,3,3);
                hrates = subplot(2,3,4);
                hrates_position = get(hrates,'Position');
                hrates_position(3) = hrates_position(3)*2;
            end
            set(hrates,'Position',hrates_position);
            h = errorbar(rates+spont_rate,stdev_of_mean, 'k','linestyle','none', 'marker','o',  'markerfacecolor','k');
            reset(gca)
            % Set post-plot axes properties
            h_axes = get(h,'parent');
            srf_hfig(2) = h_axes;
            set(h_axes,'xtick',[4 8 12 16 20 24 28 32]);%works for both 24 and 32
            ylabel(y_axis_label)
            xlabel('Location')
            %Show spontaneous spiking rate
            hspont = line([.5 length(rates)],[spont_rate spont_rate], 'linestyle','--','color','k');
            %Show threshold
            hthreshold = line([.5 length(rates)],[threshold threshold]+spont_rate, 'linestyle','--', 'color',[1 0 0]);
            legend([hspont hthreshold],{'Spontaneous Firing Rate','Threshold'});legend('boxoff')
        end
        
    end
elseif length(rates) == 5%plot the behavior related SRF
    
    locations = [1 2 3 7 8];
    locations_ticks = {'3', '2', '1', '8', '7'};
    locations_index = [3 2 1 5 4];
    [~, best_location_index] = max(rates);
    best_location = locations(best_location_index);
    number_of_peaks = NaN;
    tuning_area = NaN;
    tuning_vector = (rates*speaker_locations([1 2 3 7 8],3:5))/sum(rates-min(rates));
    rate_interpolation = [];
    
    if plot_on
        srf_hfig = figure('name',['SRF, unit ' num2str(unit_number)], 'color','white');
        h = errorbar(rates(locations_index)+spont_rate,stdev_of_mean(locations_index));
        % Set post-plot axes properties
        h_axes = get(h,'parent');
        set(h_axes,'XTick',[1:5], 'XTickLabel',locations_ticks);
        ylabel(y_axis_label)
        xlabel('Location')
        line([1 length(rates)],[spont_rate spont_rate], 'linestyle','--');
    end
    
else %not 24 or 32 or 5 locations - not an SRF file.
    disp('not an SRF file')
end

%% outputs
if length(rates) == speaker_number
datafile.best_location = best_location;
datafile.number_of_peaks = number_of_peaks;
datafile.tuning_area = tuning_area;
datafile.tuning_vector = tuning_vector;
datafile.tuning_vector_magnitude = tuning_vector_magnitude;
datafile.lateral_centroid = lateral_centroid;
datafile.best_lateral_angle = best_lateral_angle;
datafile.yz_centroid = yz_centroid;
datafile.rate_interpolation = rate_interpolation;
datafile.srf_hfig = srf_hfig;
datafile.driven_rates = rates;
end
end
