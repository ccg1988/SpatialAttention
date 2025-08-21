% compute the SRF, shift it to the center and overlap all neurons together
% combine the other 10 neurons from "SRF_inhibtion_single"
% by CCG @ 2020-05-29

clear; clc; 
% close all
load ( 'firing_rates_SRF.mat' )
% load('aa.mat')
% Unit_rates = [Unit_rates; Unit_rates_sub];
% Unit_rates_spont = [Unit_rates_spont; Unit_rates_spont_sub];
%%
N = length(Unit_rates_spont);
num_effect_2D=nan(N, 24);
num_peak=nan(N, 1);
best_location=nan(N, 1);
inhibit_ratio=nan(N, 24);
inhibit_ratio_all=nan(N, 1);
num_inhibit=nan(N, 1);
change_elevation = 0 ;

if change_elevation ==1
    Unit_rates(57, 6) = Unit_rates(57, 6) + 1;     %keep it static strongly inhibited! 
    Unit_rates(75, 4) = Unit_rates(75, 4) + 2;     %keep it static
    Unit_rates(128, 6) = Unit_rates(128, 6) + 3; %keep it static
    Unit_rates(133, 6) = Unit_rates(133, 6) + 4; %keep it static
    Unit_rates(134, 5) = Unit_rates(134, 5) + 1; %keep it static
    Unit_rates(419, 7) = Unit_rates(419, 7) + 15;
    Unit_rates(620, 8) = Unit_rates(620, 8) + 3;
    Unit_rates(620, 24) = Unit_rates(620, 24) + 3; % avoid the bottom zero line
end

for n = 1 : N
% for n = 1 : 1
    % 3 spikes sound 1 spike spont:    (3-1)/(3+1)=0.5       1 spike sound 3 spikes spont: (1-3)/(1+3)=-0.5
    % no spikes during spontaneous: (7.8-0)/(7.8+0)=1      no spikes during sound: (0-3.58)/(0+3.58)=-1
    inhibit_ratio_all(n)=sum(Unit_rates(n, :)) - 24*Unit_rates_spont(n); %44<0    0==0    605>0      
    inhibit_ratio(n, :)=(Unit_rates(n, :) - Unit_rates_spont(n)) ./ (Unit_rates(n, :) + Unit_rates_spont(n));
    temp_spon = Unit_rates(n, :) - Unit_rates_spont(n);
    [~, best_location(n)] = max (temp_spon);
    num_inhibit(n) = length(find (temp_spon<0));
    num_effect_2D(n, :) =  sign(temp_spon);
    temp_thres = Unit_rates(n, :) - max ( Unit_rates(n, :) )/2; %those speakers that pass half threshold
    num_peak(n) = length(find (temp_thres>=0));
end

if change_elevation ==1
%     ID = find ( (Unit_rates_spont >=1) & (sum(inhibit_ratio,2)<=-1) & ( num_peak<=15 ) &...
%     (best_location~=16)& (best_location~=24) ) ;
    ID = find ( (Unit_rates_spont >=1) & (sum(inhibit_ratio,2)<= 0) & (best_location~=16) & (best_location~=24) ) ;
else
    ID = find ( (Unit_rates_spont >=1) & (sum(inhibit_ratio,2)<= 0 ) ) ; %%%%%%%%%% Keep parameter %%%%%%%%%
end
num_ID= length (ID);

%%
% show the SRF of single neurons, original and modulated
% n = 653 ;
% rates = Unit_rates(n, :);
% spont_rate = Unit_rates_spont(n);
% 'T' is transform  SRF;
% 1 is only SRF, 2 is both SRF+Spikes
% [rate_interpolation] = analyze_srf_simple( rates,  spont_rate, n, num_inhibit(n), 'T', 2 );  
%%
rates_raw = zeros (num_ID, 24);
empty_id_2D = zeros (num_ID, 1);
empty_line = nan (1, 72);
for i = 1 : num_ID
    temp_id=ID(i);
    rates25 = [Unit_rates(temp_id,:)   Unit_rates_spont(temp_id)];
    max25 = max (rates25);
    % normalize the firing rates among 24 locations + spontaneous rate
    rates_raw (i, :) = Unit_rates(temp_id,:)/max25 - Unit_rates_spont(temp_id)/max25 ;
end

load('SRF.mat'); 
speaker_locations = speakers; clear speakers
spatial_resolution = 5 ;
longitude_limit = 180-spatial_resolution/2;
latitude_limit = 90-spatial_resolution/2;
nlongitude_cells = 360/spatial_resolution;
nlatitude_cells = 180/spatial_resolution;
C_x=nlongitude_cells/2; C_y=nlatitude_cells/2;

longitude = repmat(-1*longitude_limit:spatial_resolution:longitude_limit,nlatitude_cells,1);
latitude = repmat(-1*latitude_limit:spatial_resolution:latitude_limit,nlongitude_cells,1)';

distances = zeros(24, nlatitude_cells, nlongitude_cells);
for i = 1 : 24 %faster to do it this way instead of element by element. Takes .01 - .3 seconds
     distances(i,:,:) = reshape(distance(speaker_locations(i,[2 1]), [latitude(:) longitude(:)]),nlatitude_cells,nlongitude_cells);
end
 
rate_interpolation_3D = nan(num_ID, nlatitude_cells,nlongitude_cells);
for n = 1 : num_ID
% for n = 56 : 56  
    rate_interpolation = zeros(nlatitude_cells,nlongitude_cells);
    rates = squeeze(rates_raw(n,:)) ;
    for i = 1:nlatitude_cells
        for j = 1:nlongitude_cells
            [dist, index] = sort(distances(:,i,j));
            weights = (1./dist(1:2)).^2./(sum(1./dist(1:2).^2));
            rate_interpolation(i,j) = sum(weights'.*rates(index(1:2)));
        end
    end   
%     rate_interpolation = flipud(rate_interpolation);
%     rate_interpolation=rate_interpolation/max(rate_interpolation,[],'all');
    [~,peak_x] = max( (max(rate_interpolation,[],1)) );
    [~,peak_y] = max( (max(rate_interpolation,[],2)) );
    T_x=C_x-peak_x;
    T_y=C_y-peak_y;
    temp_x = circshift (rate_interpolation ,T_x, 2 ); %2 for X-demension
      
    % method 1: circular shift of Y/elevation axis
%     temp_interpolation = circshift (temp_x ,T_y, 1 ); % no empty area, rotate bottom-up 
    % method 2: do not change Y/elevation axis
    temp_interpolation=imtranslate(temp_x,[0, 0]);
    % method 3: change elevation and pat with NaN
    if change_elevation ==1
        temp_interpolation = imtranslate ( temp_x, [0, T_y] ); % empty area will be zero (==spon)
        empty_id = find(sum(temp_interpolation,2)==0);
        
        empty_sign = 0;
    if max(empty_id)==36
        empty_sign=-1;  %empty in the top
        for  e = min(empty_id) : max(empty_id)
%             temp_interpolation ( e, : ) = temp_interpolation ( min(empty_id) - 1, : );
            temp_interpolation ( e, : ) = empty_line;
        end
    elseif min(empty_id)==1 
        empty_sign=1;
        for  e = 1 : max(empty_id)
%             temp_interpolation ( e, : ) = temp_interpolation ( max(empty_id) + 1, : );
            temp_interpolation ( e, : ) = empty_line;
        end
    end
    empty_id_2D(n) = empty_sign * length(empty_id);
    
    end
    
    rate_interpolation_3D(n, :, :)=temp_interpolation;
    
end
%%
% show the interpolated plot (flipud)
if change_elevation ==1
    % rate_interpolation_3D_m=squeeze(nanmean(rate_interpolation_3D)); %1 and NaN to 1
    rate_interpolation_3D_m=squeeze(nansum(rate_interpolation_3D))/num_ID; %1 and NaN to 0.5   
else
    rate_interpolation_3D_m=squeeze( mean(rate_interpolation_3D) ); %W/O nan in this matrix
end
% normalize the peak-values to 1
 % rate_interpolation_3D_m=rate_interpolation_3D_m/max(rate_interpolation_3D_m,[],'all');

filter_sigma = 1 ;
% [4 1] anisotropic Gaussian smoothing kernels (different STD along row and column) 
rate_interpolation_smooth = imgaussfilt (rate_interpolation_3D_m , filter_sigma ); %2D isotropic Gaussian filter

pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position',[X_size*0.05 Y_size*0.2 X_size*0.45 Y_size*0.45]);
pos1 = [0.05 0.06 0.55 0.78]; % length: X-axis = Y-axis*2
subplot('Position', pos1)
% imagesc( flipud(rate_interpolation_3D_m) );
imagesc( flipud (rate_interpolation_smooth) ); hold on
[M,~] = contour( flipud (rate_interpolation_smooth), [0 0]); %0 is spontaneous
M (:,1) = [] ;
% show the contour using spontaneous firing rates
plot( M(1, :), M(2, :), '--', 'Color', rgb('white'), 'LineWidth', 1)
figure
[M_new,~] = contour( (rate_interpolation_smooth), [0 0] ); %0 is spontaneous  'ShowText','on'
M_new (:,1) = [] ; close

% M_new = flipud(M);
% M_new = rot90(M, 3);
% figure;plot( M_new(1, :), M_new(2, :), '--', 'Color', rgb('black'), 'LineWidth', 1)

xticks ([1 10 19 28 36.5 45 54 63 72])
xticklabels ( { '-180', '-135', '-90', '-45', '0', '45', '90', '135', '180' } )
yticks ([1 9 18 27 36 ])
yticklabels ( { '90', '45', '0', '45', '90' } )
title( [ 'Number: ', num2str(num_ID), '  \sigma = ',  num2str(filter_sigma) ] )
% axis off
colormap(jet)
% colorbar('Ylim',[-1 1],'location','westoutside');
bar_min = min(rate_interpolation_smooth, [] , 'all' );
bar_max = max(rate_interpolation_smooth, [] , 'all' );
colorbar('Ylim', [bar_min bar_max] ,'location','westoutside');

% show the Azimuth & Elevation value curves across peak-point
[ ~, x_id ]=max( (max(rate_interpolation_smooth,[],1)) ); %azimuth direction
% [ ~, y_id ]=max( (max(rate_interpolation_smooth,[],2)) ); %the same as following
% y_id start from 36/bottom, so y_id 14 is displayed as 23
[ ~, y_id ]=max( rate_interpolation_smooth(:, x_id) ); %elevation direction
pos2 = [0.65 0.06 0.3 0.7];
subplot('Position', pos2)
x1= 1.5 : 72.5 ;
h1 = plot(x1, rate_interpolation_smooth( y_id, : ) ); %line of azimuth     1~72   Left to Right
hold on
x2= 18+( 1.5 : 36.5 )+(18-y_id) ;
h2 = plot( x2, rate_interpolation_smooth( :, x_id ) ); %line of elevation   1~36  Down to Up
hold on
h3 = plot( x1, zeros(72,1), '--', 'color', 'k');
% xticks ([1 19 36 54 72])
% xticklabels ( { '-180', '-90', '0', '90', '180' } )
xticks ([1.5 10 19 28 36.5 45 54 63 71.5])
xticklabels ( { '-180', '-135', '-90', '-45', '0', '45', '90', '135', '180' } )
xlim([1.5 71.5]); 
y_min = min(rate_interpolation_smooth,[],'all');
y_min = round (y_min,2);
ylim([bar_min bar_max])
legend( [h1 h2 h3] , {'Azimuth', 'Elevation', 'Spontaneous'} );

tuning_mean = rate_interpolation_smooth( y_id, : );
posi_id = find(tuning_mean>0); neg_id = find(tuning_mean<0);
posi_num = num2str(length(posi_id)); 
neg_num = num2str(length(neg_id));
posi_area = round(sum(tuning_mean(posi_id),2), 2);
neg_area = round(sum(tuning_mean(neg_id),2), 2);
title( {[ '# units: ', num2str(length(ID)), ' # pos vs neg: ',  posi_num, ' ', neg_num], ...
    [' area pos vs neg: ',  num2str(posi_area), ' ', num2str(neg_area)]} )
%%
% show the SRF map based on interpolated firing rates
figure_title = 'CCG';
srf_hfig = figure('name',figure_title,'position',[30 250 600 325]);
ax = axesm('MapProjection', 'mollweid');%fournier is used by Evan's paper
setm(ax,'Grid','on','Frame','off','MlabelParallel',10,'LabelFormat','none');
axis off
spatialR = georasterref('RasterSize', size(rate_interpolation_smooth), 'Latlim', [-90 90], 'Lonlim', [-180 180]); 
geoshow(rate_interpolation_smooth, spatialR,'DisplayType', 'texturemap'); 
colormap(jet)
hemisphere_line = [-90*ones(181,1) (-90:90)'; 90*ones(181,1) (90:-1:-90)'];       
geoshow(hemisphere_line(:,2), hemisphere_line(:,1),'color',[.5 .5 .5],'linewidth',1);

M_scale =  M_new ;
M_scale(1,:) = M_scale (1,:)*spatial_resolution-(180+spatial_resolution/2);
M_scale(2,:) = M_scale (2,:)*spatial_resolution-(90+spatial_resolution/2);
M_scale(1,M_scale(1,:) == max(max(longitude))) = 180;
M_scale(1,M_scale(1,:) == min(min(longitude))) = -180;
M_scale(2,M_scale(2,:) == max(max(latitude))) = 90;
M_scale(2,M_scale(2,:) == min(min(latitude))) = -90;
% M_scale = rot90 (M_scale, 3) ;
% M_scale = flipud (M_scale, 3) ;
geoshow( M_scale(2, :), M_scale(1, :), 'color',  rgb('white') ,'linewidth',1)

colorbar('Ylim',[bar_min bar_max],'location','westoutside');
x_position = min(get(gca,'xlim'))-0.5;%shift "y_axis_label" to the left
y_position = min(get(gca,'ylim'));
text(x_position,y_position,'Spikes/s',  'HorizontalAlignment','Left', 'VerticalAlignment','Top'); 
%%
% show the 3D surf of "rate_interpolation_smooth"
figure; 
[ X , Y ] = meshgrid ( 1:72 , 1:36 ) ; %not necessary to use this
surf( X , Y , rate_interpolation_smooth, 'EdgeColor', 'flat')
% colormap(jet); 
% colorbar('Ylim',[bar_min bar_max],'location','westoutside');
xlim([1 72]); xticks ([1 19 36.5 54 72])
xticklabels ( { '-180\circ',  '-90\circ',  '0\circ',  '90\circ',  '180\circ' } )
ylim([1 36]); yticks ([1 9 18 27 36 ])
yticklabels ( { '90\circ', '45\circ', '0\circ', '45\circ', '90\circ' } )