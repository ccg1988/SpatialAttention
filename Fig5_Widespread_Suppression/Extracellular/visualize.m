clear; clc
% datafile_rates=analyze_rates('M9X',[2157 4 ],0,[1,1]); % prefer only front and back; main figure
% datafile_rates=analyze_rates('M9X',[305 4 ],0,[1,1]); % inhibition everywhere, put in supplementary figure
datafile_rates=analyze_rates('M71V',[1209 5 ],0,[1,1]); % only one speaker but larger variation, behavior session; main figure
% datafile_rates=analyze_rates('M3T',[710 5 ],0,[1,1]); % wrong
% datafile_rates=analyze_rates('M3T',[816 1 ],0,[1,1]); % Evan shown in thesis but not CC paper; main figure


%2nd value is for plotting: 1(SRF), 2(SRF+PSTH), 3(SRF+PSTH+Raster)
%3rd value is direction 1(facing front)/2(facing left)/3(facing right)
%4th value is name for matching "YYW" or "ER" or any others
[datafile_SRF, best_location, number_of_peaks, tuning_area, ~,...
    rate_interpolation, srf_hfig] = analyze_srf_beta(datafile_rates, 2 ,1 , 'c'); %'T' for tranform the SRF