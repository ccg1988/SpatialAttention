function [mV] = raw2mV(raw)

% modify the wave (mV)
% program from Jen's folder

amp_gain = 10;
system_bitdepth = 16;
scale_factor = (((10/2^(system_bitdepth-1)))/amp_gain).*10^3; %mV
mV = raw*scale_factor;