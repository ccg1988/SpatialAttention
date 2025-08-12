function r = raw_postsynaptic_activity(D1,X1,Y1,G_sum1, Th,convolution_mode)
% eq=1or9; G_sum1/2 is 29*29, D/X/Y is 5*5*29*29*29*29
sum1=zeros(1,length(X1)); 
switch  convolution_mode
    case 'Gaussian_weights'
           % for x = 1:57
           %      for y = 1:57
           %          Large_RF(x,y) = exp((-(x-29)^2)/(2*sgm^2)).*exp((-(y-29)^2)/(2*sgm^2));  
           %      end
           % end 
           %  for n = 1:size(G_sum1,2)
           %      for m = 1:size(G_sum1,2)
           %          conv_G_sum1(n,m) = sum(sum(Large_RF((30-n):(58-n),(30-m):(58-m)).* G_sum1));
           %          conv_G_sum2(n,m) = sum(sum(Large_RF((30-n):(58-n),(30-m):(58-m)).* G_sum2));
           %      end
           %  end
           %   for k = 1:length(X1)
           %        sum1(k) = conv_G_sum1(X1(k),Y1(k)); sum2(k) = conv_G_sum2(X2(k),Y2(k));    
           %   end
    case 'point_weights'    
           for k = 1:length(X1) % 5*5*29*29*29*29
                 sum1(k) = G_sum1(X1(k),Y1(k)); 
           end
    otherwise
           disp('please enter a correct convolution_mode');
           return;   
end     
      membrane_potential = D1.*sum1 - Th; % equation 9
      membrane_potential(membrane_potential<0) = 0;
      r = membrane_potential;
return






















% %%%%%%%%%%  Gaussian weight convolution, slow computation 1 %%%%%%%%%%%%
% sgm = 2;
%    for x = 1:33
%         for y = 1:33
%             Large_RF(x,y) = exp((-(x-17)^2)/(2*sgm^2)).*exp((-(y-17)^2)/(2*sgm^2));  
%             
%         end;
%     end;
%  
%     sum1=zeros(1,length(X1));
%     sum2=zeros(1,length(X1));
%     for k = 1:length(X1)  
%         sum1(k) = sum(sum(Large_RF((18-X1(k)):(34-X1(k)),(18-Y1(k)):(34-Y1(k))).* G_sum1));
%         sum2(k) = sum(sum(Large_RF((18-X2(k)):(34-X2(k)),(18-Y2(k)):(34-Y2(k))).* G_sum2));
% 
%     end



   %%%%%%%%%%  Gaussian weight convolution, slow computation 2 %%%%%%%%%%%%
% sgm = 2;
%    for x = 1:57
%         for y = 1:57
%             Large_RF(x,y) = exp((-(x-29)^2)/(2*sgm^2)).*exp((-(y-29)^2)/(2*sgm^2));  
%             
%         end;
%     end;
%  
%     sum1=zeros(1,length(X1));
%     sum2=zeros(1,length(X1));
%     for k = 1:length(X1)  
%         sum1(k) = sum(sum(Large_RF((30-X1(k)):(58-X1(k)),(30-Y1(k)):(58-Y1(k))).* G_sum1));
%         sum2(k) = sum(sum(Large_RF((30-X2(k)):(58-X2(k)),(30-Y2(k)):(58-Y2(k))).* G_sum2));
%    
%     %end




%%%%%%%%%%%  Gaussian weight convolution, fast computation  %%%%%%%%%%%%
%  sgm = 2;
%  
%    for x = 1:57
%         for y = 1:57
%             Large_RF(x,y) = exp((-(x-29)^2)/(2*sgm^2)).*exp((-(y-29)^2)/(2*sgm^2));  
%             
%         end;
%     end; 
%     
% 
%     for n = 1:size(G_sum1,2)
%         for m = 1:size(G_sum1,2)
% 
%             conv_G_sum1(n,m) = sum(sum(Large_RF((30-n):(58-n),(30-m):(58-m)).* G_sum1));
%             conv_G_sum2(n,m) = sum(sum(Large_RF((30-n):(58-n),(30-m):(58-m)).* G_sum2));
%     
%         end
%     end
% 
%      sum1=zeros(1,length(X1));
%      sum2=zeros(1,length(X1));
%    for k = 1:length(X1)
%        
%         sum1(k) = conv_G_sum1(X1(k),Y1(k));
%         sum2(k) = conv_G_sum2(X2(k),Y2(k));
%       
%    end
