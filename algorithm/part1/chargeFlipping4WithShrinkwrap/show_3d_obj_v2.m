% Used to plot 3D object
% This code is used to show 3d density.
% First created by Yun Zhao on Sep 6, 2014
% Last editted by Yun Zhao on Feb 13, 2015

% Compared with V1, this shows certain percentage of voxels.

prompts = 'Please give the obj file: \n';
obj_file = input(prompts,'s');
if isempty(obj_file)
    obj_file = 'support_hio_large.mat';
end
obj = importdata(obj_file);
% obj = abs(obj_k);
% obj = real(obj_k);
% obj = importdata('unit_cell_1A.mat');
% obj = fftshift(obj);
% obj = modulus_centered;
[s_a,s_b,s_c]=size(obj);

prompts = 'Please set the solvent fraction level (in percentage): \n';
fraction_level_percent = input(prompts);
fraction_level = fraction_level_percent/100.0;
obj_1d = reshape(obj,[1,s_a*s_b*s_c]);
obj_1d_sort = sort(obj_1d);
I_contour = obj_1d_sort(int64(fraction_level*s_a*s_b*s_c));


% I_1D_array = reshape(obj,1,s_a*s_b*s_c); % Squeeze 3d array to 1d array
% I_mean = mean(I_1D_array);
% I_sigma = std(I_1D_array);
% I_contour = sigma_level*I_sigma;

figure(5)
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
