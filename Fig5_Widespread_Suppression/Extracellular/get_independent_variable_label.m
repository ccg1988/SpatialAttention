function [ var_label ] = get_independent_variable_label(var_ind, analysis_code)

%--------------------------------------------------------------------------
% Get Name of Varying Parameter(s) (Adapted from D. Bendor)
%--------------------------------------------------------------------------

if var_ind==0
    label = 'stim #';
elseif var_ind == 1
    label = 'speaker #';
elseif var_ind == 2
    label = 'atten (dB)';
else
    switch analysis_code
        case 1 %bandwidth
            label{1}='frequency (Hz)';
        case 2  %rate-level
            label{1}='frequency (Hz)';
        case 10 %tones
            label{1}='frequency (Hz)';
            label{3}='duration (ms)';
        case 32  %clicks
            label{1}='ICI (ms)';
        case 33  %gaussian noise click
            label{1}='ICI (ms)';
        case 34  %click_test1
            label{1}='sigma';
            label{2}='ICI (ms)';
        case 35 %click_test2
            label{1}='sigma';
            label{2}='ICI (ms)';
        case 36 %rclick
            %label{1}='ICI (ms)';
            label{1}='Click Width (ms)';
            label{2}='ICI (ms)';
        case 37 %gclick
            %label{1}='sigma';
            %label{2}='ICI (ms)';
            label{1}='sigma';
            label{2}='ICI (ms)';
            label{3}='Fc (Hz)';
        case 38 %gnclick
            %label{1}='ICI (ms)';
            label{1}='sigma';
            label{2}='ICI (ms)';
        case 43  % lFM
            label{1}='CF (kHz)';
            label{2}='depth (Hz)';
        case 51 %sAM
            label{1}='center frequency (kHz)';
            label{2}='mod. depth';
            label{3}='modulation freq. (Hz)';
        case 52  % sFM
            label{1}='CF (kHz)';
            label{2}='depth (Hz)';
            label{3}='modulation freq. (Hz)';
        case 53  % sMM
            label{1}='CF (kHz)';
            label{2}='FM depth (Hz)';
            label{3}='AM depth';
            label{4}='modulation freq. (Hz)';
            label{5}='Phase diff (deg)';
        case 54 %BP noise AM
            label{1}='center frequency (kHz)';
            label{2}='half-BW (Hz)';
            label{3}='AM depth';
            label{4}='modulation freq. (Hz)';
            label{5}='randn seed';
            label{6}='randn seed';
        case 55 %WB - AM
            label{1}='AM depth';
            label{2}='modulation freq. (Hz)';
            label{3}='randn seed';
            label{4}='randn seed';
        case 60 % Noise
            label{1}='CF (kHz)';
            label{2}='bandwidth (oct.)';
            label{3}='randn seed';
            label{4}='randn seed';
        case 61 % WB-noise
            label{1}='CF (kHz)';
            label{2}='bandwidth (oct.)';
            label{3}='randn seed';
            label{4}='randn seed';
        case 62 %BP-noise
            label{1}='center frequency (kHz)';
            label{2}='bandwidth (oct.)';
        case 63 % BR-noise
            label{1}='CF (kHz)';
            label{2}='bandwidth (oct.)';
            label{3}='randn seed';
            label{4}='randn seed';
        case 64 % LP-noise
            label{1}='Cutoff freq. (kHz)';
            label{2}='bandwidth (oct.)';
            label{3}='randn seed';
            label{4}='randn seed';
        case 65 % HP-noise
            label{1}='Cutoff freq. (kHz)';
            label{2}='bandwidth (oct.)';
            label{3}='randn seed';
            label{4}='randn seed';
        case 1301 %linear_fm
            label{1}='CF (kHz)';
            label{2}='atten (dB)';
            label{3}='bandwidth';
            label{4}='duration';
            label{5}='F_start (kHz)';
            label{6}='F_stop (kHz)';
            label{7}='nrep';
            label{8}='trajectory';
            label{9}='esn';
            label{10}='mxamp';
        case 20 %two_tone2
            label{1}='freq 1 (kHz)';
            label{2}='atten 1 (dB)';
            label{3}='freq 2 (kHz)';
            label{4}='atten 2 (dB)';
            label{5}='tone 2 bw';
        case 1205
            label{1} = 'Fc(Hz)';
            label{2} = 'fRange(oct)';
            label{3} = 'Steps_per_oct';
            label{4} = 'tRange(ms)';
            label{5} = 'Bin_size(ms)';
            label{6} = 'Pip len (ms)';
            label{7} = 'Pip F(Hz)';
            label{8} = 'Pip start (ms)';
            label{9} = 'Pip BW (oct)';
            for j = 10:40
                label{j} = 'two_pip';
            end
        otherwise
            if analysis_code>=100 && analysis_code<=196 %user
                label{1} = 'stim #';
            elseif analysis_code>=600 && analysis_code<=699  %Dennis's stimuli (Comodulation)
                label{1} = 'stim #';
            elseif analysis_code>=700 && analysis_code<=799  %Dennis's stimuli (RSS)
                label{1}='SR';
                label{2}='Ramp';
                label{3}='Set#';
                label{4}='MWBEStim#';
                label{5}='Freq (Hz)';
                label{6}='BW (oct)';
                label{7}='bins/oct';
                label{8}='tones/bin';
                label{9}='Level Contrast (dB stdev)';
                label{10}='Atten Adj (dB)';
                label{11}='AM freq';
                label{12}='AM depth';
                label{13}='AM coherence';
                label{14}='FM freq';
                label{15}='FM depth';
                label{16}='FM coherence';
                label{17}='ID';
            elseif analysis_code>1100 && analysis_code<1200 %pblaster
                if analysis_code<1111
                    label{1}='delay (ms)';
                elseif analysis_code<1151
                    label{1}='ICI (ms)';
                elseif analysis_code<1161
                    label{1}='tone density (Hz)';
                elseif analysis_code<1171
                    label{1}='2nd tone freq (Hz)';
                elseif analysis_code>1170
                    label{1}='harmonic complex code';
                end
                
                if analysis_code<1111
                    label{2}='iterations';
                elseif analysis_code<1113
                    label{2}='sigma (ms)';
                elseif analysis_code<1123
                    label{2}='decay (ms)';
                elseif analysis_code<1133
                    label{2}='width (ms)';
                elseif analysis_code<1143
                    label{2}='slope (oct/10 ms)';
                elseif analysis_code<1161
                    label{2}='bandwidth (Hz)';
                elseif analysis_code>1170
                    label{2}='# of harmonics';
                end
                
                label{3}='jitter (%)';
                
                if analysis_code<1111
                    label{4}='noise cf (Hz)';
                elseif analysis_code<1111
                    label{4}='center frequency (Hz)';
                elseif analysis_code<1151
                    label{4}='frequency';
                elseif analysis_code<1161
                    label{4}='center frequency (Hz)';
                elseif analysis_code>1170
                    label{4}='fundamental frequency (Hz)';
                end
                
                if analysis_code<1111
                    label{5}='noise bw (Hz)';
                elseif analysis_code == 1141
                    label{5} = 'intensity_slope (dB/width)';
                    label{6} = 'width (ms)';
                elseif analysis_code>1170 && analysis_code<1180
                    label{5}='stimulus attenuation (dB)';
                    label{6}='atten norm';
                    label{7}='freq1 string' ;
                    label{8}='freq2 string' ;
                    label{9}='jitter string';
                    label{10}='harmonics range';
                elseif analysis_code>1180 && analysis_code<1190
                    label{5}='stimulus attenuation (dB)';
                    label{6}='noise center frequency (Hz)';
                    label{7}='noise bandwidth (Hz)';
                    label{8}='S/N ratio (dB)';
                    label{9}='noise atten (dB)';
                    label{10}='atten norm';
                    label{11}='freq1 string' ;
                    label{12}='freq2 string' ;
                    label{13}='jitter string';
                    label{14}='harmonics range';
                end %if
            elseif analysis_code>1700 && analysis_code<1800
                label{1} = 'freq (kHz)';
                label{2} = 'bandwidth (oct)';
                label{3} = 'preduration (ms)';
                label{4} = 'dur/cyc';
                label{5} = 'on1 (ms)';
                label{6} = 'off1 (ms)';
                label{7} = 'att1';
                label{8} = 'onset1';
                label{9} = 'offset1';
                label{10} = 'on2 (ms)';
                label{11} = 'off2 (ms)';
                label{12} = 'att2';
                label{13} = 'onset2';
                label{14} = 'offset2';
            elseif analysis_code == 2101 %CIstim
                label{7} = 'current (uAmp)'; %7
                label{6} = 'pulse rate (hz)'; %6
                label{10} = 'electrode #'; %10
                label{11} = 'pulse width (usec)'; %11
                label{12} = 'interphase gap (usec)'; %12
                label = label{var_ind};
                %                 label = label(var_ind-3)
            end %if
            %                         label = label{var_ind-2};
    end %switch
    
end %if
var_label = label;

end
