% Change the spherical degree (azi/ele/r) to the Cartesian degree (XYZ)
% Input matrix is 24*3 or 32*3   1st is 1:N, 2nd is elevation, 3rd is azimuth   
% Out matrix is 24*5 or 32*5     1st is azimuth, 2nd is elevation, 3rd-4th-5th is Y-X-Z
% tuning vector is [left(negative)/right(positive) front(positive)/rear(negative) up(positive)/down(negative)]
clear;clc
load('speakers1_32_YYW.mat');
azi=deg2rad(speakers(:,2));
ele=deg2rad(speakers(:,3));
value=ones(size(speakers,1),1);%=1, because all speakers are 1m distance
[y,x,z]=sph2cart(azi,ele,value);  %sqrt(x^2+y^2+z^2)=1

x(abs(x)<0.001)=0;y(abs(y)<0.001)=0;z(abs(z)<0.001)=0;
new_SRF=[speakers(:,2) speakers(:,3) x y z];
% isequal(new_SRF,speakers)
%%
%This is used to practice the transformation of YXZ and Azi/Ele 
clear; clc
load('SRF.mat')
[azimuth,elevation,r]=cart2sph(speakers(:,4),speakers(:,3),speakers(:,5)); %notice Y-X-Z, similar as sph2cart
% This is the sentence from "analyze_srf_beta.m"
% [THETA,PHI,~] = cart2sph(tuning_vector(2),tuning_vector(1),tuning_vector(3));
azimuth=rad2deg(azimuth);%-180deg to 180deg
azimuth_full=azimuth;
azimuth_full(azimuth<0)=azimuth(azimuth<0)+360;
elevation=rad2deg(elevation);
%%
% clear; clc
datafile=analyze_rates('M9X',[2157 4 ],0,[1,1]);
% load('SRF.mat')%speaker information
% in case of 32 speakers
rates = datafile.rates; rates_new = rates(1:24); datafile.rates = rates_new;
rates_std = datafile.stdev_of_mean; rates_std_new=rates_std(1:24); datafile.stdev_of_mean = rates_std_new;
%2nd value is for plotting: 1(SRF), 2(SRF+PSTH), 3(SRF+PSTH+Raster)
%3rd value is direction 1(facing front)/2(facing left)/3(facing right)
%4th value is name for matching "YYW" or "ER" or any others
[datafile, best_location, number_of_peaks, tuning_area, ~,...
    rate_interpolation, srf_hfig] = analyze_srf_beta(datafile, 2 ,1 , 'c'); %'T' for tranform the SRF
% azimuth=rad2deg(azimuth); azimuth_full=azimuth; azimuth_full(azimuth<0)=azimuth(azimuth<0)+360;
% elevation=rad2deg(elevation); 
% p=datafile.p_values;p_sig01=length(find(p<0.01));p_sig001=length(find(p<0.001));
% rate_interpolation = flipud(rate_interpolation);
%%
spont_rate=0;
[tuning_area, tuning_vector, tuning_vector_magnitude, azimuth, elevation] = analyze_srf_alpha(rates,spont_rate,speakers);
%%
% shift the SRF
% peak_x=max(rate_interpolation, [],'all');
C_x=36; C_y=18;
[~,peak_x] = max( (max(rate_interpolation,[],1)) );
[~,peak_y] = max( (max(rate_interpolation,[],2)) );
T_x=C_x-peak_x;
T_y=C_y-peak_y;
temp_x = circshift (rate_interpolation ,T_x, 2 ); %2 for X-demension
% temp = circshift (temp_x , T_y , 1 ); %1 for Y-demension; should not shift
temp = imtranslate ( temp_x, [0, T_y] );
figure;imagesc(temp); colormap(jet)
% figure;imagesc(rate_interpolation); colormap(jet)
%%
clear; clc
datafile=analyze_rates('M71V',[524 4],0,[1,1]);
load('SRF.mat')%speaker information
[tuning_vector, tv_magnitude, azimuth, elevation] = analyze_srf_alpha(datafile.rates,datafile.spont_rate,speakers);
%Display the firing rates of each repeat of two/three selected speakers
%Used for the dynamical SRF changes across time
peak1=5;%Back  Ele 0
peak2=1;%Front Ele 67 ofr NO9
% peak3=1;
figure;
% rates_diff=abs(datafile.rep_rates(:,peak2)-datafile.rep_rates(:,peak1));
% ax3=plot(rates_diff,'k--o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on
ax1=plot(datafile.rep_rates(:,peak1),'k--*','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');%[0.5,0.5,0.5]
hold on;
ax2=plot(datafile.rep_rates(:,peak2),'k--s','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
ylim([0 50])
% hold on
% ax3=plot(datafile.rep_rates(:,peak3),'k--o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
% ylabel('Spikes/Second');lgd=legend([ax1 ax2 ax3],{num2str(peak1),num2str(peak2),'diff'});title(lgd,'Speaker Location')
ylabel('Spikes/Second');lgd=legend([ax1 ax2],{num2str(peak1),num2str(peak2)});title(lgd,'Speaker Location')
% ylabel('Spikes/Second');lgd=legend([ax1 ax2 ax3],{num2str(peak1),num2str(peak2),num2str(peak3)});title(lgd,'Speaker Location')
R_duration=datafile.stim*(datafile.pre_stim+datafile.post_stim)+sum(datafile.stim_dur);%ms
xlabel(['#Repeat:     ',num2str(R_duration/1000),' Second/Repeat']); 
title(['File: ',datafile.datafile_num,'  Channel: ', num2str(datafile.ch)])
%%
%Calculate the time gap between two files
%Used for the dynamical SRF changes across time
data_raw1st=M9X0677; 
data_raw2nd=M9X0680;
str1st = data_raw1st.datetime;% '28-Mar-2012 11:10:46'
str2nd = data_raw2nd.datetime;% '28-Mar-2012 11:10:46'
t1st= datevec(str1st,'dd-mmm-yyyy HH:MM:SS');
t2nd= datevec(str2nd,'dd-mmm-yyyy HH:MM:SS');
t_absolute=etime(t2nd,t1st);%time-diff of 2 files' created time
t_duration=(data_raw1st.data(end,4)+data_raw1st.post_stimulus_record_time*1000)/1000000;%stimulus duration in first file: us to s
t_relative=t_absolute-round(t_duration);%time between the end of 1st file and the start of 2nd file
