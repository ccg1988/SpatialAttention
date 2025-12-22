function [ tuning_area, tuning_vector_magnitude] = analyze_srf_half_YW( rates, spon, side, latlim, plot_on )

if size(rates,1)>size(rates,2)
   rates = rates';
end

figure_title = ( 'CCG' ) ;

if strcmp(side, 'Right') == 1
%     srf = load('SRF_half_R.mat');
    srf = load('Speakers_no_4_right.mat');
elseif strcmp(side, 'Left') == 1
%     srf = load('SRF_half_L.mat');
    srf = load('Speakers_no_4_left.mat');
end

if length(rates) == 14
    srf = load('Speakers_no_4_11_left.mat');
end

speaker_locations = srf.speakers;
speaker_number=size(speaker_locations,1);
spatial_resolution = 5;
longitude_limit = 180-spatial_resolution/2;
latitude_limit = 90-spatial_resolution/2;
nlongitude_cells = 360/spatial_resolution;
nlatitude_cells = 180/spatial_resolution;

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
    
    [rate_max, best_location] = max(rates);
    
%     threshold = .5*(rate_max+spon)-spon;
    threshold = .5*(rate_max);
    
    
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
    
    tuning_area = sum(interp_cell_areas(rate_interpolation > threshold));
    
    
    %*********************************Tuning vector and angles*********************************
    %tuning_vector = mean([rates' rates' rates'].*speaker_locations(:,3:5))/mean(rates);---Wrong
    % old tuning vector is biased by speaker distribution (analyze_SRF, L: 231)
    % tuning_vector = (rates*speaker_locations(:,3:5))/sum(rates-min(rates));
    
    interp_vectors = zeros(nlatitude_cells*nlongitude_cells,3);
    [interp_vectors(:,2), interp_vectors(:,1), interp_vectors(:,3)] = ...
        sph2cart(pi*longitude(:)/180,pi*latitude(:)/180, ...
        (rate_interpolation(:)+spon).*interp_cell_areas(:));
    tuning_vector = 1.066*sum(interp_vectors)/sum((rate_interpolation(:)+spon).*interp_cell_areas(:));
    %1.066 is an adjustment so that a response to only one location has a magnitude of 1;
    [~,~,tuning_vector_magnitude] = cart2sph(tuning_vector(2),tuning_vector(1),tuning_vector(3));
    tuning_vector_magnitude = min(tuning_vector_magnitude,1);
    
%**************************************************************Plot SRF*****==1/2/3******************************************
    if plot_on  
        if plot_on == 1
            srf_hfig = figure('name',figure_title,'position',[30 250 600 425]);
        elseif plot_on == 2
            srf_hfig = figure('name',figure_title,'position',[30 250 525 525]);
        elseif plot_on == 3
            srf_hfig = figure('name',figure_title,'position',[30 250 820 525]);
        end
        
%         rate_interpolation_raw=rate_interpolation;
        rate_interpolation (1:18, :) = [] ; %remove the lower part of SRF
        latlim_1 = [0 90] ;
%         latlim_2 = [-30 90] ;
        latlim_2 = latlim ;
        spatialR = georasterref('RasterSize', size(rate_interpolation), 'Latlim', latlim_1 , 'Lonlim', [-180 180]); 
        %spatialref.GeoRasterReference object
        ax = axesm('MapProjection', 'mollweid', 'MapLatLimit',latlim_2);%fournier is used by Evan's paper
        %Show the angles of azimuth & elevation lines
%         setm(ax,'Grid','on','Frame','off','MeridianLabel','On','MlabelParallel',10, ...
%             'ParallelLabel','On','LabelFormat','none');
        setm(ax,'Grid','on','Frame','off','MlabelParallel',10,'LabelFormat','none');
        axis off
        
        if plot_on == 3
            set(gca,'position',[.02 .5 .64 .5]);
        elseif plot_on == 2
            set(gca,'position',[.02 .5 1 .5]);
        else
            set(gca,'position',[0.02 0 0.97 1]);
        end
        
        %Show the inter-polated spike rates
        geoshow(rate_interpolation, spatialR,'DisplayType', 'texturemap'); 
        colormap(jet)
        absolute_maximum = max(abs(rates));
        if length(unique(rates)) > 1 %Avoiding error
            set(gca,'Clim',[-1*absolute_maximum absolute_maximum]);
        end
        
        %Show the white-circle hemisphere line
        hemisphere_line = [-90*ones(181,1) (-90:90)'; 90*ones(181,1) (90:-1:-90)'];     
        geoshow(hemisphere_line(:,2), hemisphere_line(:,1),'color',[.5 .5 .5],'linewidth',1);
        
        location_numbers = 1:speaker_number;
        locations_index = zeros(1,speaker_number);
        
        passive_locations = ones(1,speaker_number);
        number_of_passive = speaker_number;    
        
        locations = struct('Geometry','Point', ...
            'lat',mat2cell(speaker_locations(~~passive_locations,2),ones(1,number_of_passive)), ...
            'long',mat2cell(speaker_locations(~~passive_locations,1), ones(1,number_of_passive)), ...
            'Name',cellstr(num2str(location_numbers(~~passive_locations)')));
        
        %Show speaker locations
        locations (17:end)=[]; %do not show the lower 8 speakers
        geoshow(locations,'Marker','o', 'MarkerFaceColor',[.5 .5 .5], 'MarkerEdgeColor','k', 'MarkerSize',7);
             
        %%
        if plot_on == 1
%             text(.4, -1.0, {['Lateral Centroid: ' num2str(lateral_centroid,2)], ['Vertical Centroid: ' num2str(yz_centroid,3)]});
            text(.4, -1.05, {['Tuning Area: ' num2str(tuning_area,2)], ['Vector Magnitude: ' num2str(tuning_vector_magnitude,2)]});
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
            x_position = min(get(gca,'xlim')) - 0.6 ; %shift "y_axis_label" to the left
            y_position = min(get(gca,'ylim'));
            y_axis_label = 'Spikes/S';
            text(x_position,y_position,y_axis_label,  'HorizontalAlignment','Left', 'VerticalAlignment','Top');       
        end       
               
    end
    
end
