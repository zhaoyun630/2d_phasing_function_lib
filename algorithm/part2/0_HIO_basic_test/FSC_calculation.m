% This script is used to calculate Fourier shell correlation (FSC)

% Created by Yun Zhao on May 12, 2016

% obj1 = importdata('3rdx.mat');
obj1 = importdata('inputID2_newObj.mat');
obj2 = importdata('runID_5_obj.mat');

% obj2 = obj1;
bin_num = 10;
res_max = 2;
FCC = zeros(1,bin_num+1);

[s_a, s_b, s_c] = size(obj1);

F1 = fftshift(fftn(obj1));
F2 = fftshift(fftn(obj2));

F1F2 = real(F1.*conj(F2));
% F1F2 = abs(F1).*abs(F2);
F1F2_sq1 = abs(F1).^2;
F1F2_sq2 = abs(F2).^2;

mask = zeros(s_a,s_b,s_c);
% The following matrix should be calculated for only once
for i=1:s_a
    for j=1:s_b
        for k=1:s_c
            mask(i,j,k)=ceil(bin_num*norm([(2*i-s_a-1)/(s_a-1),(2*j-s_b-1)/(s_b-1),(2*k-s_c-1)/(s_c-1)]));
        end
    end
end
            
% mask = ceil(bin_num*mask);


for i=1:bin_num+1
    mask_i = (mask == i-1) ;
    a1 = F1F2.*mask_i;
    a2 = F1F2_sq1.*mask_i;
    a3 = F1F2_sq2.*mask_i;
    FCC(i) = sum(a1(:))/(sqrt(sum(a2(:)))*sqrt(sum(a3(:))));
end

x=linspace(0,res_max,bin_num+1)/bin_num;
plot(x,FCC);
axis([0 0.2 0 1.1])
