function [lateral_centroid, best_lateral_angle, yz_centroid] = horizontal_pole_transform(rates, plot_on)

%This function bases tuning calculations based on an alternate coordinate
%system where longitudes are parallel and latitudes converge at the
%contralateral and ipsilateral most points on the sphere. In this case,
%longitude will be -90 to +90 (centered at the front), and latitude will
%vary from -180 to 180, with 0 at the front. For the projection, map will be
%cut and unwrapped from 180 (rear).

%% Setup spatial grid
srf = load('SRF.mat');
speaker_locations = srf.speakers;
spatial_resolution = 5;

longitude_limit = 180-spatial_resolution/2; %#ok<*NASGU>
latitude_limit = 90-spatial_resolution/2;
nlongitude_cells = 360/spatial_resolution;
nlatitude_cells  = 180/spatial_resolution;

%% Remap speaker locations
x_component = speaker_locations(:,4);
y_component = speaker_locations(:,5);
z_component = speaker_locations(:,3);

[theta, phi, rho] = cart2sph(x_component, y_component, z_component);
theta = theta*180/pi;
phi = -1*phi*180/pi;

speaker_locations_transformed = [phi theta];

%% Calculate rate interpolation
longitude_list = -1*longitude_limit:spatial_resolution:longitude_limit;
latitude_list = -1*latitude_limit:spatial_resolution:latitude_limit;
longitude = repmat(longitude_list,nlatitude_cells,1);
latitude = repmat(latitude_list,nlongitude_cells,1)';

distances = zeros(24, nlatitude_cells, nlongitude_cells);

for i = 1:24
    distances(i,:,:) = reshape(distance(speaker_locations_transformed(i,[1 2]), ...
        [latitude(:) longitude(:)]),nlatitude_cells,nlongitude_cells);
end

% Testing
% rates = zeros(1,24);
% rates(2) = 1;
% plot_on = 2;
% spont_rate = 0;

rate_interpolation = zeros(nlatitude_cells,nlongitude_cells);


for i = 1:nlatitude_cells
    for j = 1:nlongitude_cells
        [dist index] = sort(distances(:,i,j));
        weights = (1./dist(1:2)).^2./(sum(1./dist(1:2).^2));
        rate_interpolation(i,j) = sum(weights'.*rates(index(1:2)));
    end
end

mean_lateral_rates = [fliplr(mean(rate_interpolation,2)) ...
    ((-90+spatial_resolution/2):spatial_resolution:(90-spatial_resolution/2))']; %change "wrev" to "fliplr" by CCG

%% This algorithm can get a spont subracted rate, or non-spont subtracted
%*This may seem a little confusing. Centroid calculations are sensitive to
%the dynamic range and minimum firing rate when using a non-circular axis,
%or treating a circular axis like it isn't one. I'm trying to reduce that
%effect, but it's a little arbitrary.
mean_lateral_rates_inhibition = mean_lateral_rates(:,1);
mean_lateral_rates_inhibition(mean_lateral_rates(:,1) > 0) = 0;
suppression_adjusted_mean_lateral_rates = ...
    fliplr(-1*mean_lateral_rates_inhibition) + mean_lateral_rates(:,1);   %change "wrev" to "fliplr" by CCG
suppression_adjusted_mean_lateral_rates(mean_lateral_rates(:,1) < 0) = 0;

lateral_centroid = sum(suppression_adjusted_mean_lateral_rates.*mean_lateral_rates(:,2))/ ...
    sum(suppression_adjusted_mean_lateral_rates);
[~, best_lateral_angle_index] = max(mean_lateral_rates(:,1));
best_lateral_angle = spatial_resolution*best_lateral_angle_index - spatial_resolution/2 - 90;
%% Calculate yz centroid
%Latitude weights so averaging is done accoring to area of cell
latitude_weight = zeros(1,nlatitude_cells);
for i = 1:nlatitude_cells
    latitude_up = latitude(i,:) + spatial_resolution/2;
    latitude_down = latitude(i,:) - spatial_resolution/2;
    longitude_left = longitude(i,:) + spatial_resolution/2;
    longitude_right = longitude(i,:) - spatial_resolution/2;
    latitude_weight(i) = sum(areaquad(latitude_down,longitude_right,latitude_up,longitude_left));
end
latitude_weight = repmat(latitude_weight',1,nlongitude_cells);
mean_yz_rates = sum(rate_interpolation.*latitude_weight,1)';

% Calculate vectors for each angle. Use vectors so that the centroid
% calculation doesn't get biased towards the center of whatever range is
% used for the degree value (like in the lateral centroid calculation).
[yz_vectors(:,1) yz_vectors(:,2)] = pol2cart((pi*longitude_list/180)',mean_yz_rates);

%get mean vector, convert to degrees
yz_vector_sum = sum(yz_vectors);
yz_centroid = 180*cart2pol(yz_vector_sum(:,1), yz_vector_sum(:,2))/pi;

mean_yz_rates = [mean_yz_rates longitude_list'];
%%
if plot_on
    figure;
    spatialR = georasterref('RasterSize', size(rate_interpolation), ...
        'Latlim', [-90 90], ...
        'Lonlim', [-180 180]); %% spatialref.GeoRasterReference object
    ax = axesm('MapProjection', 'fournier');
    %         setm(ax,'Grid','on','Frame','off','MeridianLabel','On','MlabelParallel',10, ...
    %             'ParallelLabel','On','LabelFormat','none');
    setm(ax,'Grid','on', ...
        'Frame','off', ...
        'MlabelParallel',10, ...
        'LabelFormat','none');
    axis off
    %srf_hfig(2) = gca; not sure why this is even output
    if plot_on == 3
        set(gca,'position',[.02 .5 .64 .5]);
    elseif plot_on == 2
        set(gca,'position',[.02 .5 1 .5]);
    else
        set(gca,'position',[0.02 0 1 1]);
    end
    
    geoshow(rate_interpolation, spatialR, ...
        'DisplayType', 'texturemap');
    absolute_maximum = max(abs(rates));
    if length(unique(rates)) > 1 %Avoiding error
        set(gca,'Clim',[-1*absolute_maximum absolute_maximum]);
        
        hbar = colorbar('Ylim',[min(rates) max(rates)], ...
            'location','westoutside');
    end
    
    location_numbers = 1:24;
    passive_locations = ones(1,24);
    number_of_passive = 24;
    
    locations = struct('Geometry','Point', ...
        'lat',mat2cell(speaker_locations_transformed(~~passive_locations,1),ones(1,number_of_passive)), ...
        'long',mat2cell(speaker_locations_transformed(~~passive_locations,2), ones(1,number_of_passive)), ...
        'Name',cellstr(num2str(location_numbers(~~passive_locations)')));
    geoshow(locations,'Marker','o', ...
        'MarkerFaceColor',[0 0 0], ...
        'MarkerEdgeColor','k', ...
        'MarkerSize',7);
    textm([locations(:).lat]*.8,[locations(:).long],{locations(:).Name},'FontSize',12)
    if plot_on == 2
        %%
        subplot(5,1,4)
        plot(mean_lateral_rates(:,2),mean_lateral_rates(:,1))
        line([lateral_centroid lateral_centroid], ...
            [min(mean_lateral_rates(:,1)) max(mean_lateral_rates(:,1))])
        set(gca,'XTick', [-90 -45 0 45 90]);
        legend('x rate','Location','EastOutside')
        subplot(5,1,5)
        plot(mean_yz_rates(:,2),mean_yz_rates(:,1))
        line([yz_centroid yz_centroid], ...
            [min(mean_yz_rates(:,1)) max(mean_yz_rates(:,1))])
        set(gca,'XTick', [-180 -135 -90 -45 0 45 90 135 180]);
        legend('yz rate', 'Location','EastOutside')
    end
end
end
