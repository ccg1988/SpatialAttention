%exponent takes 1,2,3,4,5. tomo03/18/09

%save the vectorized population parameters, except for the variable
%thresholds, due to the memory limitation. tomo3/17/09


%clear all
  
% current_time=fix(clock);
% disp(sprintf('start at %d:%d on %d/%d/%d',current_time(4),current_time(5),current_time(1),current_time(2),current_time(3)))

D1_array=[1.0];
X1_array=1:29;
Y1_array=1:29;

N=length(D1_array)*length(X1_array)*length(Y1_array);

ni = 1; %neuron index
D1 = zeros(1,N,'single');
X1 = zeros(1,N,'single');
Y1 = zeros(1,N,'single');

for l=1:length(Y1_array)
    for m=1:length(X1_array)
            for o=1:length(D1_array)
         
                  Y1(ni)=Y1_array(l);
                  X1(ni)=X1_array(m);
                  D1(ni)=D1_array(o);
                  ni = ni + 1;

            end
        end
    end

% disp('vect_params created')
% 
% 
% save ('Vector_params.mat','alpha','beta','visX0','visY0','vestX0','vestY0','expo','alpha_array','beta_array','visX0_array',...
%     'visY0_array','vestX0_array','vestY0_array','expo_array');
% disp(sprintf('end at %d:%d on %d/%d/%d',current_time(4),current_time(5),current_time(1),current_time(2),current_time(3)))
% 
% 
% 
% 
