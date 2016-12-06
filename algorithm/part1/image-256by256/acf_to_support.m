% This piece of code is used to get a support for autocorrelation from the
% autocorrelation calculation of support.

% acf_of_support = importdata('acf_obj_1x.mat');

acf = importdata('acf_half.mat');

% [size_a,size_b]= size(acf_of_support);
% acf_1d = reshape(acf_of_support,[1,size_a*size_b]);
I_max = max(acf(:));
I_mean = mean(acf(:));
I_cut = 0.0000001*I_max;
mask = acf > I_cut;
support_acf = mask.*ones(size(acf));



% for i = 1:size_a
%     for j = 1:size_b
%         if acf_of_support(i,j)>=I_cut
%             support_for_acf(i,j)=1;
%         end
%     end
% end
% 
% save('support_from_acf_2x.mat','support_acf');
save('support_from_half_acf_1x.mat','support_acf');

imshow(support_acf)