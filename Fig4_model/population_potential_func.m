% call>>> CreateVectorParams, StimlusVector2Matrix, 
% call>>> raw_presymaptic_activity(eq=567), saturation_func(eq=8), raw_postsynaptic_activity(eq=1or9)
function population_potential = population_potential_func(stim_D1, stim_X1, stim_Y1, sgm)
CreateVectorParams; %returns post-synaptic neuron's parameters: D1,X1,Y1 & D2,X2,Y2
saturation_func = @hypoEXPO_f;
expo = 0.5; % the number used only in the case for the hypoEXPO presynaptic saturation function   
post_synaptic_convolution_mode = 'point_weights' ; %'Gaussian_weights'
Th = 0;% a number used for the post-synaptic threshold             
stimulus_matrix1= [stim_D1 stim_X1 stim_Y1 sgm]; % 2by3 matrix, [c1,15, 15; c2, 15, 15] auditory
G_sum1=saturation_func(raw_presymaptic_activity(stimulus_matrix1), expo); % 29*29 auditory>>>>>>>same as final outputs
resp  = raw_postsynaptic_activity(D1,X1,Y1,G_sum1, Th,post_synaptic_convolution_mode); % 1*29*29
population_potential(1:length(D1_array), 1:length(X1_array), 1:length(Y1_array))...
    =reshape(resp,length(D1_array), length(X1_array), length(Y1_array));
population_potential = squeeze(population_potential);
        
       