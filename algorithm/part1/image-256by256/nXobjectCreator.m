% create a n-X bigger obj
% Written by Yun in 3/13/2016

obj = importdata('obj_cd_256.mat');
support = importdata('support_cd_256.mat');
[s1, s2]=size(obj);

% size of new obj
oversampling_ratio = 2; % oversampling in each dimension
n1 = oversampling_ratio*s1;
n2 = oversampling_ratio*s2;

obj_nX = zeros(n1,n2);
support_nX = zeros(n1,n2);
obj_nX((n1-s1)/2+1:(n1+s1)/2,(n2-s2)/2+1:(n2+s2)/2)= obj;
support_nX((n1-s1)/2+1:(n1+s1)/2,(n2-s2)/2+1:(n2+s2)/2)= support;

objfileName = strcat('obj_cd_',int2str(oversampling_ratio),'x.mat');
supportfileName = strcat('support_cd_',int2str(oversampling_ratio),'x.mat');
save(objfileName,'obj_nX');
save(supportfileName,'support_nX');
figure(1)
imshow(obj_nX)
figure(2)
imshow(support_nX)
