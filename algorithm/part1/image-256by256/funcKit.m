% function kit
% created by Yun Zhao 3/13/2016
% "cd" in file name means "cat and dog"

% Get support from image
% supportFile = 'support_cd_256.jpg';
% im1 = im2double(importdata(supportFile));
% support = im1(:,:,1);
% save('support_cd_256.mat','support');

%%%%
% Create support based on size information
cd_2x = importdata('obj_cd_2x.mat');
[n1,n2] = size(cd_2x);
s1 = n1/2;
s2 = n2/2;
support_2x_square = zeros(n1,n2);
support_2x_square((n1-s1)/2+1:(n1+s1)/2,(n2-s2)/2+1:(n2+s2)/2)=1;
save('support_2x_square.mat','support_2x_square');

%%%
% Get half of image (for acf demonstration);
acf = importdata('Patt_fun.mat');
[n1,n2]=size(acf);
s1 = n1/2;
s2 = n2/2;
acf_half = acf((n1-s1)/2+1:(n1+s1)/2,(n2-s2)/2+1:(n2+s2)/2);
save('acf_half.mat','acf_half');
imagesc(acf_half);
axis equal
axis off

% Calculate Patterson function for crystal.
obj = importdata('obj_cd_256.mat');
F_complex = fft2(obj);
F_modules = abs(F_complex);
I = F_modules.^2;
patterson_function_cd_256 = fftshift(abs(ifft2(I)));
save('patterson_function_cd_256.mat','patterson_function_cd_256');
imagesc(patterson_function_cd_256);
axis equal
axis off

% Show a range of values 
obj1 = importdata('difpat_512.mat');
objMax = max(obj1(:));
upperValue = 0.1*objMax;
bottomValue = 0.00001*objMax;
mask1 = obj1 < upperValue;
mask2 = obj1 > bottomValue;

obj1 = mask1.*obj1 + (1-mask1)*upperValue;
obj1 = mask2.*obj1 + (1-mask2)*bottomValue;
save('diff_pattern_512_demo.mat','obj1');
imagesc(obj1)
axis equal
axis off

iptsetpref('ImshowBorder','tight');
set(gca,'position',[0 0 1 1],'units','normalized')
saveas(gca,'test.jpg');
% imwrite(obj1,'test.tif');


