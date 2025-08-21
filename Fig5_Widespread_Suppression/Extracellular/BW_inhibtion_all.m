% Compute the side-band inhibition of spectral tuning curves
% Similar as the "SRF_inhibtion_all"
% by CCG @ 2020-05-30

clear; clc; close all
load ( 'firing_rates_BW.mat' )
total_N = 649;
center = 21 ; % 1:20 21 22:41
S = 41 ;
shift_tuning_2D = nan (649, S) ;
inhibit_ratio = nan (649, S) ;
rates_max = nan (649, 1);
for n = 1 : total_N
% for n = 18 : 19
    stimuli = Unit_stimuli{n} ;
    if length(stimuli) == S
        raw_spike = Unit_spike_raw{n} ;
        spike_spont = Unit_spontaneous(n);
        spike_mean = mean (raw_spike, 2) ;
        [ BF_peak_value, BF_peak_id ] = max (spike_mean) ;
        shift_step = center - BF_peak_id ;
        % not use the peak of two ends
        if ( BF_peak_id~=1 && BF_peak_id~=S )
           shift_tuning = circshift ( spike_mean ,shift_step );
           tuning42 = [ shift_tuning; spike_spont ];  max42 = max(tuning42); rates_max(n) = max42;
           temp = ( shift_tuning - spike_spont ) / max42;
           shift_tuning_2D ( n, : ) = temp ;          
           inhibit_ratio ( n, : ) = ( spike_mean - spike_spont ) ./ ( spike_mean + spike_spont );
           
        end
    end
end
%%
% show the tuning curve of single neuron
pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position',[X_size*0.05 Y_size*0.4 X_size*0.5 Y_size*0.35], 'Visible','on');

N = 24 ;

subplot(1,2, 1)
raw_spike = Unit_spike_raw{N};
% h1 = plot(raw_spike); hold on
h2 = plot(mean(raw_spike, 2),  '--*', 'LineWidth', 1, 'Color', 'k'); hold on
h3 = plot( ones(S, 1)*Unit_spontaneous(N) ); box off
legend( [h2 h3] , {'Tuning mean', 'Spont'} );
xlim([1 S])
xticks ([ 1 11 21 31 41 ])
xticklabels ( { '2', '4', '8', '16', '32' } )
xlabel('Pure tone freq. (kHz)'); ylabel('Firing rate (Spikes/S)')
title( [ 'Example unit: ', num2str(N), '  spont: ',  num2str(Unit_spontaneous(N)), ' spikes/s' ] )

subplot(1,2, 2)
h1 = plot( shift_tuning_2D (N,:), '--*', 'LineWidth', 1, 'Color', 'k'); box off
% hold on
% h2 = plot( ones(S, 1)*Unit_spontaneous(N)/rates_max(N) );
legend( {'Tuning curve'} );
xlim([1 S])
xticks ([ 1 11 21 31 41 ])
xticklabels ( { '-2', '-1', '0', '1', '2' } )
xlabel('Distance from best freq. (Octave)'); ylabel('Normalized firing rate (Spikes/S)')
title( [ 'Example unit: ', num2str(N), '  spont: ',  num2str(Unit_spontaneous(N)), ' spikes/s' ] )
%%
figure;
ID = find ( ( Unit_spontaneous>=1 ) & ( sum(inhibit_ratio,2)<= 0 ) ) ;
ID_tuning = shift_tuning_2D (ID, :);
% tuning_std = std(ID_tuning,[],2)/length(ID);
tuning_std = std(ID_tuning,[],2) ;
tuning_mean = mean ( ID_tuning, 1 );
% tuning_mean = imgaussfilt (tuning_mean,1); %smooth tuning curve
% tuning_mean = smooth(tuning_mean,0.05); %smooth tuning curve
% h1 = errorbar ( tuning_mean,tuning_std );
% shadedErrorBar( 1:S, tuning_mean, tuning_std,'g');
h1 = plot(tuning_mean);
hold on
h2 = plot( ones(S, 1)*0, '--' ); box off
legend( [h1 h2] , {'Tuning curve', 'Spontaneous'} );
xlim([1 S])
xticks ([ 1 11 21 31 41 ])
xticklabels ( { '-2', '-1', '0', '1', '2' } )

posi_id = find(tuning_mean>0); neg_id = find(tuning_mean<0);
posi_num = num2str(length(posi_id)); 
neg_num = num2str(length(neg_id));
posi_area = round(sum(tuning_mean(posi_id),2), 2);
neg_area = round(sum(tuning_mean(neg_id),2), 2);
title( [ '# units: ', num2str(length(ID)), ' # pos vs neg: ',  posi_num, ' ', neg_num, ...
    ' area pos vs neg: ',  num2str(posi_area), ' ', num2str(neg_area)] )