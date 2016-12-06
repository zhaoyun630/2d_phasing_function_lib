
% Used to import image from jpg file
% Written by Yun Zhao 3/13/2016


%location of cat
c1 = 30;
c2 = 30;

% location of duck
d1 = 128 + 10;
d2 = 128 + 10;

% input image
% duck = 'duck-logo.jpg';
% cat = 'kit-kat.jpg';
duck = 'duck1.jpeg';
cat = 'cat4.jpg';

% import duck image
im1 = im2double(importdata(duck));
im2 =abs(im1(:,:,1)-1);
im3 = imresize(im2,0.4);
[a1,a2] = size(im3);

%import cat image
im4 = im2double(importdata(cat));
im5 = im4(:,:,1);
x1 = (im5>0.84);
im5 = (1-x1).*im5;
im6 = imresize(im5,0.1);
[b1,b2] = size(im6);

obj = zeros(256,256);
obj(d1:d1+a1-1,d2:d2+a2-1) = im3;
obj(c1:c1+b1-1,c2:c2+b2-1) = im6;

save('obj_cd_256.mat','obj');
imshow(obj)
imwrite(obj,'obj_cd_256.jpg');
% % Cat 2
% im4 = im2double(importdata(cat));
% im5 = im4(50:380,:,1);
% x1 = (im5>0.999);
% im5 = (1-x1).*im5;
% im6 = imresize(im5,[50,50]);
% [a2,b2] = size(im6);

% % Cat 3
% im4 = im2double(importdata(cat));
% im5 = im4(150:850,30:540,1);
% x1 = (im5>0.99);
% im5 = (1-x1).*im5;
% im6 = imresize(im5,0.2);
% [a2,b2] = size(im6);