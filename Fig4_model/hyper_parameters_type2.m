clear;clc;
close all
width_field = 29 ; %entire field of image
N = width_field*width_field;
% rng("default") using generator 'twister' and seed 0
seed = 0 ; rng(seed, "twister"); 
% 6+3 generators: https://www.mathworks.com/help/matlab/ref/rng.html
norm_sgm_HP = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20] ;
norm_drift_peak_HP = [1, 2, 3, 4] ; % the maximum drift
attn_drift_peak_HP = [1, 2, 3, 4, 6, 8, 10, 12] ; % the maximum drift

attn_smaller_all = [] ;
ratio_all = [] ;
nHigh_all = [] ;
nMid_all = [] ;
normSGM_all = [] ;
norm_all = [] ;
attn_all = [] ;

for normSGM_HP = 1 : length(norm_sgm_HP)
    for norm_HP = 1 : length(norm_drift_peak_HP)
        for attn_HP = 1 : length(attn_drift_peak_HP)

stim_amp_drift = 8.75 + [1,2,3,4,5] ; % 8.75; fixed; control firing rate within 100spk/s
norm_sgm = norm_sgm_HP(normSGM_HP) ; % only one value, no need loop
norm_drift=randi(norm_drift_peak_HP(norm_HP), 5,1)-1; % peak_drift, #num
attn_drift=rand(50,1)*attn_drift_peak_HP(attn_HP); % attn_drift=round(attn_drift,2);

div_sgm_all = rand( length(stim_amp_drift) * length(norm_drift) * length(attn_drift), 1 );
peak_activity_in = nan(length(stim_amp_drift), length(norm_drift), length(attn_drift));
peak_activity_out = nan(length(stim_amp_drift), length(norm_drift), length(attn_drift));
num = 0 ; % total number of neurons
% norm_sgm is determined by hyper-para. but is constant in loop for each session
for s= 1:length(stim_amp_drift) % determine each point in each session
for n= 1:length(norm_drift)
for a= 1:length(attn_drift)

stim_D1 = 2^stim_amp_drift(s) ;
stim_X1 = 15 ;
stim_Y1 = 15 ;
stim_sgm = 2 ;
stim_PP = population_potential_func(stim_D1, stim_X1, stim_Y1, stim_sgm);
stim_div_expo = 2;   
stim_PP = stim_PP.^stim_div_expo;

attn_D1 = 1 ;
attn_X1 = stim_X1 - attn_drift(a);
attn_Y1 = stim_Y1 - attn_drift(a) ;
attn_sgm = 2 ;
attn_PP = population_potential_func(attn_D1, attn_X1, attn_Y1, attn_sgm);
attn_div_expo = 2;   
attn_PP = attn_PP.^attn_div_expo;

norm_D1 = 512 ;
norm_X1 = stim_X1 - norm_drift(n) ;
norm_Y1 = stim_Y1 - norm_drift(n) ;
norm_PP = population_potential_func(norm_D1, norm_X1, norm_Y1, norm_sgm);
norm_div_expo = 2;  
norm_PP = norm_PP.^norm_div_expo;

num = num + 1 ;
div_sgm = 32*div_sgm_all(num) ;
stim_attn2norm_in = (stim_PP*attn_PP)./(div_sgm^norm_div_expo + norm_PP); 
stim_attn2norm_out = (stim_PP*1)./(div_sgm^norm_div_expo + norm_PP); 
peak_activity_in(s, n, a)   = stim_attn2norm_in  (stim_X1,stim_Y1);
peak_activity_out(s, n, a) = stim_attn2norm_out(stim_X1,stim_Y1);

end
end
end

fr_in_raw = peak_activity_in(:); fr_out_raw = peak_activity_out(:);
fr_keep = find((fr_in_raw>1)&(fr_out_raw>1));
fr_in = fr_in_raw(fr_keep); fr_out = fr_out_raw(fr_keep);
attn_smaller = 100*length(find(fr_out>fr_in))/ length(fr_keep) ;
attn_smaller_all = [attn_smaller_all, attn_smaller];
xScaled = fr_out / max(fr_out);
nHigh = sum(xScaled >= 0.7 & xScaled <= 1.0); nHigh_all = [nHigh_all, nHigh];
nMid  = sum(xScaled >= 0.2 & xScaled <= 0.5); nMid_all = [nMid_all, nMid];
ratio = nHigh / max(nMid,eps); ratio_all = [ratio_all, ratio];
normSGM_all = [normSGM_all, norm_sgm] ;
norm_all = [norm_all, norm_drift_peak_HP(norm_HP)-1] ;
attn_all = [attn_all, attn_drift_peak_HP(attn_HP)] ;

%
if (attn_smaller>10) && (attn_smaller<30) && (nHigh>50) && (nMid>50) && (ratio>0.4)
pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
f = figure('position',[X_size*0.05 Y_size*0.1 X_size*0.9 Y_size*0.8],'visible','off'); 
subplot(1,2, 1)
scatter(fr_out, fr_in, 'x', 'black'); hold on 
plot(linspace(-5, 100), linspace(-5, 100), '--', 'Color', 'black','LineWidth', 2)

axis square;xlabel('Passive spikes/second'); ylabel('Hits spikes/second')
title(['Attention decreased=', num2str(round(attn_smaller)), '%'])
xlim([-5, 100]); ylim([-5, 100]); xticks([0 50 100]); yticks([0 20 40 60 80 100]);

subplot(1,2, 2);
scatter(xScaled, fr_in./fr_out, 'x', 'black');
xlim([0, 1]); xticks([0 0.2 0.4 0.6 0.8 1]); ylim([0, 4]); yticks([0 1 2 3 4]);
axis square; xlabel('Scaled rate passive'); ylabel('Scaled increase hits')
title(['Norm. sigma=', num2str(norm_sgm),' ratio=',num2str(round(ratio,2))])

temp_file_name = strcat('normSGM_',num2str(norm_sgm), ...
    '_norm_',num2str(norm_drift_peak_HP(norm_HP)-1),...
    '_attn_',num2str(attn_drift_peak_HP(attn_HP)), ...
    '_attnSmall_',num2str(round(attn_smaller)), ...
    '_ratio_',num2str(round(ratio,2)),'.jpg');
saveas(f, temp_file_name)

end

        end
    end
end

save('stimuli_results','seed','stim_amp_drift','norm_sgm_HP','norm_drift_peak_HP', ...
    'attn_drift_peak_HP','attn_smaller_all', 'ratio_all', 'nHigh_all', 'nMid_all', ...
    'normSGM_all', 'norm_all', 'attn_all')

% '_s_',num2str(stim_amp_drift(s)), ... % last stimulus in this session
% '_n_',num2str(norm_drift(n)), ...
% '_a_',num2str(round(attn_drift(a),1)), ...