% This piece of code is used to plot density distribution
% First created by Yun Zhao on Sep 6,2014
% Last editted by Yun Zhao on Sep 6, 2014

prompt = 'Please give the obj file: \n';
obj_file = input(prompt,'s');
if isempty(obj_file)
    obj_file = 'unit_cell_1A.mat'; % default file, if not input is given
end
obj_1 = importdata(obj_file);
[s_a1,s_b1,s_c1]=size(obj_1);

prompt = 'Is it a triple cell(y/n)? \n';
fprintf('If it is not a triple cell along c, then the program will \n');
fprintf('assume it is a single unit cell \n');
cell_type = input(prompt,'s');
if cell_type == 'y';
    obj = obj_1(:,:,s_c1/3+1:2*s_c1/3);
else
    obj = obj_1;
end


[s_a,s_b,s_c]=size(obj);
n_obj= s_a*s_b*s_c; % number of voxels in obj file.
n_sampling = 1000;

rho_1D_array = reshape(obj,1,n_obj); % Squeeze 3d array to 1d array
rho_1D_array_log = log(rho_1D_array);
rho_sigma = std(rho_1D_array);
rho_min = min(rho_1D_array);
rho_max = max(rho_1D_array);

rho_0 = (rho_max-rho_min)/n_sampling;
n_count = zeros(1,n_sampling+1);
x_I_sigma = linspace(rho_min,rho_max,n_sampling+1)/rho_sigma;
x_I_abs = linspace(rho_min,rho_max,n_sampling+1);
rho_sort = sort(rho_1D_array);
x_I_voxel = linspace(1,n_obj,n_obj);

for i = 1:n_obj
    f_ratio = (rho_1D_array(i)-rho_min)/rho_0;
    f_floor = floor(f_ratio);
    f_ceil = f_floor+1;
    f_fraction = f_ratio - f_floor;
    n_count(f_ceil) = n_count(f_ceil) + 1;
end



% to shift the graph, as F_000 is just from guess or estimation
% Here I place the peak at the origin.
% [n_count_max,n_max_location]=max(n_count);
% rho_peak_location = x_I_abs(n_max_location);

% Now save as a modified object.
% density_shift = x_I_abs(n_max_location);
% obj_new = obj-x_I_abs(n_max_location);
% save('obj_density_shift.mat','obj_new');

% now shift it to zero
% x_I_sigma = x_I_sigma - x_I_sigma(n_max_location);
% x_I_abs = x_I_abs - x_I_abs(n_max_location);
% rho_sort = rho_sort - rho_sort(n_max_location);




figure(1);
plt1=plot(x_I_sigma,n_count,'b');
% xlim([0,ceil(rho_max/rho_sigma)]);
xlabel('charge density with unit   \Delta\rho','FontSize', 15,'FontWeight','bold');
ylabel('counts','FontSize', 15,'FontWeight','bold');
title('charge density distribution within unit cell','FontSize', 15,'FontWeight','bold');
%  saveas(plt1,'density_distribution_sigma_imag_ifftshift.tif');
 % saveas(plt1,'density_distribution_sigma_real_ifftshift.tif');
% saveas(plt1,'density_distribution_sigma_real_fftshift.tif');
% saveas(plt1,'density_distribution_sigma.tif');

figure(2);
plt2=plot(x_I_abs,n_count,'r');
% xlim([0,ceil(rho_max/rho_sigma)]);
xlabel('charge density','FontSize', 15,'FontWeight','bold');
ylabel('voxel counts','FontSize', 15,'FontWeight','bold');
title('charge density distribution within unit cell','FontSize', 15,'FontWeight','bold');
%  saveas(plt2,'density_distribution_abs_imag_ifftshift.tif');
%  saveas(plt2,'density_distribution_abs_real_ifftshift.tif');
% saveas(plt2,'density_distribution_abs_real_fftshift.tif');
% saveas(plt2,'density_distribution_abs.tif');

figure(3)
plt3=plot(x_I_voxel,rho_sort,'m');
xlim([1,ceil(n_obj*1.1)]);
xlabel('pixel number','FontSize', 15,'FontWeight','bold');
ylabel('charge density','FontSize', 15,'FontWeight','bold');
title('sorted charge density over voxel number','FontSize', 15,'FontWeight','bold');
%  saveas(plt3,'density_distribution_sorted_voxel_imag_ifftshift.tif');
%   saveas(plt3,'density_distribution_sorted_voxel_real_ifftshift.tif');
% saveas(plt3,'density_distribution_sorted_voxel_real_fftshift.tif');
% saveas(plt3,'density_distribution_sorted_voxel.tif');
