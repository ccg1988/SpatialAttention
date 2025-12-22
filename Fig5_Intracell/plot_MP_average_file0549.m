clear; clc; close all
Monkey_ID = 'M117A';
Ch_Num = 6 ;
Onset_time = 0; %para. used by Sheng
Offset_time = 50; %collected spikes after sound offset
range = 37 ; %length of removed data points of either side of spike     75 is 3ms   37 is 1.5ms
File_Name = 'M117A0549'; 

File_Dir_spike = 'D:\SU\Yunyan\Spikes_Intra\';
wd = cd;
eval(['cd ' File_Dir_spike ';']);
eval(['x=' File_Name ';']);               
eval(['cd ' wd ';']);    
clear wd;
pre_stim = x.pre_stimulus_record_time ;  % in mS
post_stim = x.post_stimulus_record_time ; % in mS
stim_dur = x.stimulus_ch1 (1,5)  ;
spike_data = x.data;
n_reps = max(spike_data(:,2));
tmax = pre_stim+stim_dur+post_stim;

File_Dir_ad2 = 'D:\SU\Yunyan\ad2\';
Ufile_ad2 = strcat(File_Name, '.ad2');
wd = cd;
eval(['cd ' File_Dir_ad2 ';']);
[data_raw, sr] = get_ad2(Ufile_ad2);           
eval(['cd ' wd ';']);    
clear wd;
data = data_raw{1, 1}; clear data_raw
[~, N_stim] = size (data);
SR = sr(1); clear sr
stim_pre_time= round ( pre_stim*SR/1000 ); %1 to 4883
stim_dur_time= round ( (stim_dur+Offset_time)*SR/1000 ); %4884 to 4883+2441
[N_repeats, N_length] = size (data{1});
MP_sound = nan(n_reps,N_stim);
MP_spon = nan(n_reps,N_stim);
MP_trace_wSpike = nan(N_length, n_reps, N_stim);
MP_trace_woSpike = nan(N_length, n_reps, N_stim);
for s = 1 : N_stim
    temp_stim = data{s};
    % [N_repeats, N_length] = size (temp_stim);
    for r = 1 : N_repeats
%     for r = 1 : 1
        trace = temp_stim (r,:);
        trace = raw2mV (trace) ; % change to mV
        MP_trace_wSpike(:, r, s) = trace ;
        trials = find ((spike_data(:,2)==r)&(spike_data(:,1)==s));
        spike_data_time_raw=spike_data(trials, 4);
        spike_data_time_raw(spike_data_time_raw<1)=[];
        spike_data_time = round ( 1.025*SR*spike_data_time_raw/1000000 );
        n_spike = length (spike_data_time);
        trace_new = trace ;
        for n = 1 : n_spike
            time_point = spike_data_time(n);
            if time_point-range>0 && time_point+range<length(trace)
                trace_new(time_point-range : time_point+range) = trace_new(time_point-range)/2+trace_new(time_point+range)/2 ;
            elseif time_point-range<=0 % starting segment
                trace_new(1:time_point) = trace_new(time_point+range) ; %keep some distance from spike peak      
            elseif time_point-range>=length(trace) %ending segment
                trace_new(time_point:length(trace)) = trace_new(time_point-range) ; %keep some distance from spike peak
            end
        end
        MP_trace_woSpike(:, r, s) = trace_new ;
        MP_spon(r,s) = mean( trace_new(1 : stim_pre_time) );
        MP_sound(r,s) = sum( trace_new(stim_pre_time+1 : stim_pre_time+stim_dur_time)-MP_spon(r,s) ); %subtract baseline MP
    end
    
end
%%
Color1 = [153,112,193]/256; %violet
Color2 = [185,141,62]/256; %brown
Color3 = [100,168,96]/256; %green
Color4 = [204,84,94]/256; %pink-red
N_cut = 801  ; % visualization only
mean_woSpike_all = nan(N_length, N_stim) ;
for s = 1 : N_stim

select_id = s ;
pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position',[X_size*0.05 Y_size*0.4 X_size*0.5 Y_size*0.5], 'Visible','off');

% N_length = N_length ; %pre_stim = round(pre_stim-1000*N_cut/SR);

select_id_wSpike = squeeze(MP_trace_wSpike(:, :, select_id));
select_id_woSpike = squeeze(MP_trace_woSpike(:, :, select_id));
select_id_woSpike = select_id_woSpike - median(select_id_woSpike);

S_start = round(pre_stim/1000*SR) ;
S_end = round((pre_stim+stim_dur)/1000*SR);
patch_col= rgb('LightGray'); 
patch_x = [ S_start S_end S_end S_start ] ;
subplot(2,1,1)
allPoints = select_id_wSpike(:);
patch_y = [min(allPoints) min(allPoints) max(allPoints) max(allPoints)] ;
patch( patch_x , patch_y ,[ 1,1,1] , 'EdgeColor',patch_col,'FaceColor',patch_col); hold on
h1 = plot(select_id_wSpike(:,1),'LineWidth',1, 'color', Color1); hold on
h2 = plot(select_id_wSpike(:,2),'LineWidth',1, 'color', Color2);
h3 = plot(select_id_wSpike(:,3),'LineWidth',1, 'color', Color3); 
h4 = plot(select_id_wSpike(:,4),'LineWidth',1, 'color', Color4); 
legend ( [h1 h2 h3 h4] , {['Trial: ', num2str(1)], ['Trial: ', num2str(2)], ['Trial: ', num2str(3)], ['Trial: ', num2str(4)]}, 'location', 'eastoutside' )
xlim([N_cut N_length]); xlabel('time(ms)'); xlabel('Time (ms)'); ylabel('mV'); ylim([min(allPoints) max(allPoints)]);
xticks([ stim_pre_time S_end N_length]);
xticklabels({ num2str(pre_stim) num2str(pre_stim+stim_dur) num2str(tmax) });
title(['Stimulus: ', num2str(select_id)])

subplot(2,1,2)
allPoints = select_id_woSpike(:);
patch_y = [-5 -5 5 5] ;
patch( patch_x , patch_y ,[ 1,1,1] , 'EdgeColor',patch_col,'FaceColor',patch_col); hold on

h1 = plot(select_id_woSpike(:,1),'LineWidth',1, 'color', Color1, 'LineWidth', 0.1); hold on
h2 = plot(select_id_woSpike(:,2),'LineWidth',1, 'color', Color2, 'LineWidth', 0.1); %h2=h2-mean(h2);
h3 = plot(select_id_woSpike(:,3),'LineWidth',1, 'color', Color3, 'LineWidth', 0.1); %h3=h2-mean(h3);
h4 = plot(select_id_woSpike(:,4),'LineWidth',1, 'color', Color4, 'LineWidth', 0.1); %h4=h2-mean(h4);
h5 = plot(mean(select_id_woSpike,2),'LineWidth',1, 'color', 'k', 'LineWidth', 2); 
mean_woSpike_all(:, s) = mean(select_id_woSpike,2);
legend ( [h1 h2 h3 h4 h5] , {['Trial: ', num2str(1)], ['Trial: ', num2str(2)], ['Trial: ', num2str(3)], ['Trial: ', num2str(4)], 'Mean'}, 'location', 'eastoutside' )
xlim([N_cut N_length]); xlabel('time(ms)'); xlabel('Time (ms)'); ylabel('mV'); 
ylim([-5 5]);
xticks([ stim_pre_time S_end N_length ]);
xticklabels({ num2str(pre_stim) num2str(pre_stim+stim_dur) num2str(tmax) });

temp_file_name = strcat(File_Name, '_Stim_',num2str(select_id), '.pdf');
% saveas(gcf, temp_file_name)
exportgraphics(gcf, temp_file_name, 'ContentType','vector', 'BackgroundColor','none')

end
%%
pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position',[X_size*0.05 Y_size*0.4 X_size*0.4 Y_size*0.3], 'Visible','on');

step = 2.5 ;
cmap = colormap(parula(N_stim));

patch_col= rgb('LightGray'); 
patch_x = [ S_start S_end S_end S_start ] ;
allPoints = select_id_wSpike(:);
patch_y = [-4 -4 38 38] ;
patch( patch_x , patch_y ,[ 1,1,1] , 'EdgeColor',patch_col,'FaceColor',patch_col); hold on

for s = 1 : N_stim

select_id = s ;
plot(mean_woSpike_all(:, s)+step*(s-1), 'Color',cmap(s,:)); hold on

end
xlim([N_cut N_length]); xlabel('Time (ms)', 'FontSize', 14); ylabel('Stimulus', 'FontSize', 14); 
box off
yticks([ 0 20 35 ]); yticklabels({ num2str(1) num2str(9) num2str(15) });
xticks([ stim_pre_time S_end N_length ]);
xticklabels({ num2str(pre_stim) num2str(pre_stim+stim_dur) num2str(tmax) });
xlim([N_cut N_length-105]);

temp_file_name = strcat(File_Name, '_All_stimuli_averaged.pdf') ;
exportgraphics(gcf, temp_file_name, 'ContentType','vector', 'BackgroundColor','none')