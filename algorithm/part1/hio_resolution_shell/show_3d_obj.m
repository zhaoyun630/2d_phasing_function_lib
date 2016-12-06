% Used to plot 3D object
% This code is used to show 3d density.
% Last editted by Yun Zhao on Sep 6, 2014

% prompts = 'Please give the obj file: \n';
% obj_file = input(prompts,'s');
% obj = importdata(obj_file);
% obj = fftshift(obj);
obj = importdata('gkxy.mat');
[s_a,s_b,s_c]=size(obj);

% prompts = 'Please set the sigma level: \n';
% sigma_level = input(prompts);
sigma_level =3.5;

I_1D_array = reshape(obj,1,s_a*s_b*s_c); % Squeeze 3d array to 1d array
I_mean = mean(I_1D_array);
I_sigma = std(I_1D_array);
I_contour = sigma_level*I_sigma;

figure(4)
fv = isosurface(obj,I_contour);
p=patch(fv);
set(p,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
view(3); % axis tight
axis([1 s_a 1 s_b 1 s_c])
% axis equal
camlight 
lighting gouraud
box on
