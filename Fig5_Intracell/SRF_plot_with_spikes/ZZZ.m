% show the raster plot and spike-rate plot
% by CCG @ 2020-07-08

clear; clc
Monkey_ID = 'M117A';
Ch_Num = 6 ;
Onset_time = 0; %para. used by Sheng
Offset_time = 50; %collected spikes after sound offset

File_Name = 'M117A0549'; %should use 'M43s5522'
[~, rates_all, spon, ~] = plot_efile_rate_one_files_YW (Monkey_ID, File_Name, Ch_Num, Onset_time, Offset_time); %

% show the SRF after get "rates" and "spon" from previous cell
rates = rates_all - spon ;
side = 'Left' ; %electrode recording side
latlim = [0 90] ;
plot_on = 1 ;
[ tuning_area, tuning_vector_magnitude] = analyze_srf_half_YW ( rates, spon, side, latlim, plot_on );