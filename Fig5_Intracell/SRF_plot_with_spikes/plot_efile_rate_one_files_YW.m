function [rel_driven_spcounts, mean_driven_rates, mean_spont_rate, side] = plot_efile_rate_one_files_YW (Monkey_ID,File_Name1,...
    Ch_Num,Onset_time,Offset_time)

% 

if strcmp(Monkey_ID,'M50p') == 1
    File_Dir = 'D:\SU\Sheng\M50p\m50P_mfiles\';
elseif strcmp(Monkey_ID,'M43s') == 1
    File_Dir = 'D:\SU\Sheng\M43s\m43S_mfiles\';
elseif strcmp(Monkey_ID,'M43q') == 1
    File_Dir = 'D:\SU\Sheng\M43q\m43Q_mfiles\';
elseif strcmp(Monkey_ID,'M117A') == 1
    File_Dir = 'D:\SU\Yunyan\Spikes\';
end

wd = cd;
eval(['cd ' File_Dir ';']);
try
    eval(['x1=' File_Name1 ';']);               
catch
    disp(strcat('Can not process ',File_Name1,'.'));
    eval(['cd ' wd ';']);
    clear wd;
end
eval(['cd ' wd ';']);    
clear wd;

rate_structure1 = EFile_rate_processing(x1,Ch_Num,Onset_time,Offset_time);
side = rate_structure1.side;
spont_rates1 = rate_structure1.spont_rates;
spont_rates1 = spont_rates1.*1000;
driven_spcounts1 = rate_structure1.driven_spcounts ;
analysis_window_length1 = rate_structure1.analysis_window_length;
driven_rates1 = driven_spcounts1/analysis_window_length1*1000;

mean_spont_rate1 = mean(mean(spont_rates1));
mean_driven_rates1 = mean(driven_rates1,2);
std_driven_rates1 = std(driven_rates1,0,2);
se_driven_rates1 = std_driven_rates1/sqrt(size(driven_rates1,2));

mean_spont_rate = ( mean_spont_rate1 ) ;
mean_driven_rates = ( mean_driven_rates1)  ; %not minus spontaneous
driven_spcounts = driven_spcounts1 ; %not minus spontaneous
rel_driven_spcounts = rate_structure1.rel_driven_spcounts ;
Rep_total = rate_structure1.rep;

pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position', [X_size*0.42 Y_size*0.28 X_size*0.3 Y_size*0.3] );
if length(mean_driven_rates1)==15 
    h1=plot(1:15,mean_driven_rates1); hold on
    set(h1,'Color','k','LineWidth',2,'LineStyle','None','Marker','o','MarkerSize',8);
    errorbar(1:15,mean_driven_rates1,se_driven_rates1, 'LineStyle','none','Color','k'); hold on
    % set(h2,'Color','k','LineWidth',1.5,'LineStyle','-'); %do not show the connectted line
elseif length(mean_driven_rates1)==14 
    h1=plot(1:14,mean_driven_rates1); hold on;
    set(h1,'Color','k','LineWidth',2,'LineStyle','None','Marker','o','MarkerSize',8);
    errorbar(1:14,mean_driven_rates1,se_driven_rates1, 'LineStyle','none','Color','k'); hold on
end

legend boxoff;
plot( [1 15] , [mean_spont_rate mean_spont_rate] , 'k--' , 'LineWidth' , 2 );
h=gca;
if length(mean_driven_rates1)==15
    set(h,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15],'XTickLabel',[1 2 3 5 6 7 8 9 10 11 12 13 14 15 16],'FontName','Arial','FontSize',12);
elseif length(mean_driven_rates1)==14
    set(h,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14],'XTickLabel',[1 2 3 5 6 7 8 9 10 12 13 14 15 16],'FontName','Arial','FontSize',12);
end
set(h,'Box','off');
h_xlabel = xlabel('Speaker #');
h_ylabel = ylabel('Firing Rate (spikes/s)');
% ylim ( [0 max(mean_driven_rates1)] )
set(h_xlabel,'FontName','Arial','FontSize',14);
set(h_ylabel,'FontName','Arial','FontSize',14);
title_text = strcat(Monkey_ID,'-',num2str(x1.unit_number),  '    File: #', File_Name1(6:end), '   Total reps: ', num2str(Rep_total)  );
h_title = title ( { [title_text] [' '] } ) ;
set(h_title,'FontName','Arial','FontSize',14);
g=gcf;
set(g,'Color','White');
% legend ( [h1] ,  [ File_Name1]  )