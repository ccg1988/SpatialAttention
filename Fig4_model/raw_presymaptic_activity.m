function G_sum = raw_presymaptic_activity(stimulus_matrix)
% inputs:   [intensity, center-X, center-Y, sigma]
% outputs: 29*29
% sgm = 10 ; % sigma in the denominitor of equation 5 (default=2); larger==broader SF, 0.5=2*1+1pixels, 1=2*3+1pixels, 2=2*5+1pixels

intensity = stimulus_matrix(1);
% width_Gaussian = stimulus_matrix(2); % always 15
width_field = 29 ; % width_Gaussian*2-1; % always 29
X = stimulus_matrix(2) ;
Y = stimulus_matrix(3) ;
sgm = stimulus_matrix(4); % 
G_sum = nan(width_field, width_field);
for x = 1:width_field
    for y = 1:width_field
        G_sum(x,y) = intensity*exp((-(x-X)^2)/(2*sgm^2)).*exp((-(y-Y)^2)/(2*sgm^2));
    end
end

return;
