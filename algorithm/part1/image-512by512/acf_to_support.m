% This piece of code is used to get a support for autocorrelation from the
% autocorrelation calculation of support.

% acf_of_support = importdata('acf_obj_1x.mat');

acf_of_support = importdata('gkxy.mat');

[size_a,size_b]= size(acf_of_support);
acf_1d = reshape(acf_of_support,[1,size_a*size_b]);
I_max = max(acf_1d);
I_mean = mean(acf_1d);
I_cut = 0.1*I_max;
support_for_acf = zeros(size_a,size_b);

for i = 1:size_a
    for j = 1:size_b
        if acf_of_support(i,j)>=I_cut
            support_for_acf(i,j)=1;
        end
    end
end
% 
% save('support_for_acf_1x.mat','support_for_acf');

imshow(support_for_acf)