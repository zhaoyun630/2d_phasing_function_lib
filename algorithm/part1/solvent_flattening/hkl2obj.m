% load hkl file and do following things:
% 1, rearrage it into 3D volumn array in reciprocal space.
% 2, do inverse Fourier Transform to get obj
% 3, pad zeros around the cell to get a larger cell.

% First created by Yun Zhao in Aug 29, 2014
% Last editted by Yun Zhao in Jan 23, 2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start. Give hkl files, including phases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 prompts = ' Please give a input hkl file (with all positive and negative hkl values): \n';
 hkl_file_name = input(prompts,'s');
 M=dlmread(hkl_file_name); % for the hkl file, you may need to delete the 
                           %  first and last line and save it as txt.
 hkl = int32(M(:,1:3)); % save hkl in a new matrix and convert it into integers.
                         % convert float hkl index to integer.
 F = M(:,4);  % extract intensity
%  F_000 = 50*max(F); % Here I use the maximum F to approaximate F_000.
F_000=0;
 F1 = 100*F/norm(F);
 phase = M(:,5);
 F_comp = F.*exp(1j.*phase*pi/180);
 % B =100* abs(log(F))/(norm(abs(log(F)))); %log scale
 [L_M, W_M]=size(M);
%  
 h_max = max(hkl(:,1));
 k_max = max(hkl(:,2));
 l_max = max(hkl(:,3));
% h_max = 28;
% k_max = 28;
 
 
 Rec_lattice = zeros(2*h_max+1,2*k_max+1,2*l_max+1);
 
 % Create a mask for futhur R factor calculation
 Rec_lattice_mask = zeros(2*h_max+1,2*k_max+1,2*l_max+1);
 
 for i=1:L_M
     h = int32(hkl(i,1)+h_max+1);
     k = int32(hkl(i,2)+k_max+1);
     l = int32(hkl(i,3)+l_max+1);
     Rec_lattice(h,k,l) = F_comp(i);
     Rec_lattice_mask(h,k,l) = 1;
 end
 
 Rec_lattice(h_max+1,k_max+1,l_max+1)=F_000;
 save('rec_lattice_3A.mat','Rec_lattice');
 save('rec_lattice_mask_3A.mat','Rec_lattice_mask');


 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Do inverse Fourier transform to get obj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  Rec_lattice_shift = fftshift(Rec_lattice);
%  unit_cell = abs(ifftn(Rec_lattice_shift));
%  %unit_cell_1 = fftshift(unit_cell);
% %  unit_cell_complex = ifftn(Rec_lattice);
% %  unit_cell_complex_shift = fftshift(unit_cell_complex);
% %  unit_cell = abs(unit_cell_complex_shift);
% %  save('unit_cell.mat','unit_cell');

 unit_cell_complex = ifftn(ifftshift(Rec_lattice));%for any FFT, the first number [111] in matrix is always (000) in real space.
 unit_cell_abs = fftshift(abs(unit_cell_complex));
 unit_cell_real = fftshift(real(unit_cell_complex));
 unit_cell_imag = fftshift(imag(unit_cell_complex));
 save('unit_cell_real_ifftshift_3A.mat','unit_cell_real');
%  save('unit_cell_imag_ifftshift.mat','unit_cell_imag');
 save('unit_cell_abs_3A.mat','unit_cell_abs');

%  unit_cell_complex = ifftn(fftshift(Rec_lattice));
%  unit_cell = real(unit_cell_complex);
%  save('unit_cell_real_fftshift.mat','unit_cell');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shift the density peak at zero.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[s_a,s_b,s_c]=size(unit_cell_real);
n_obj = s_a*s_b*s_c;
rho_1d = reshape(unit_cell_abs,1,n_obj);
rho_1d_sort = sort(rho_1d);
rho_min = min(rho_1d);
rho_max = max(rho_1d);

n_sampling = 1000;
rho_0 = (rho_max-rho_min)/n_sampling;
n_count = zeros(1,n_sampling+1);
x_I_abs = linspace(rho_min,rho_max,n_sampling+1);
x_I_voxel = linspace(1,n_obj,n_obj);

for i = 1:n_obj
    f_ratio = (rho_1d(i)-rho_min)/rho_0;
    f_floor = floor(f_ratio);
    f_ceil = f_floor+1;
    f_fraction = f_ratio - f_floor;
    n_count(f_ceil) = n_count(f_ceil) + 1;
end

[n_count_max,n_max_location]=max(n_count);
rho_shift = x_I_abs(n_max_location);

unit_cell_shift = unit_cell_real - rho_shift;
x_I_abs_shift = x_I_abs - rho_shift;
rho_1d_sort_shift = rho_1d_sort - rho_shift;

% save('unit_cell_shift_3A.mat','unit_cell_shift');
% save('rho_shift_3A.mat','rho_shift');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check density distribution for shifted obj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
plt2=plot(x_I_abs_shift,n_count,'r');
% xlim([0,ceil(rho_max/rho_sigma)]);
xlabel('charge density','FontSize', 15,'FontWeight','bold');
ylabel('voxel counts','FontSize', 15,'FontWeight','bold');
title('charge density distribution within unit cell','FontSize', 15,'FontWeight','bold');
%  saveas(plt2,'density_distribution_real_ifftshift_3A.tif');
% saveas(plt2,'density_distribution_real_fftshift.tif');
saveas(plt2,'density_distribution_abs.tif');

figure(3)
plt3=plot(x_I_voxel,rho_1d_sort_shift,'m');
xlim([1,ceil(n_obj*1.1)]);
xlabel('pixel number','FontSize', 15,'FontWeight','bold');
ylabel('charge density','FontSize', 15,'FontWeight','bold');
title('sorted charge density over voxel number','FontSize', 15,'FontWeight','bold');
 saveas(plt3,'density_distribution_sorted_voxel_abs_ifftshift_3A.tif');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  % Now create a triple cell
%  [s_a,s_b,s_c] = size(unit_cell);
%  
%  % Triple in a,b,c direction
%  obj_triple_abc= zeros(3*[s_a,s_b,s_c]);
%  support_abc = zeros(3*[s_a,s_b,s_c]);
%  obj_triple_abc((s_a+1):2*s_a,(s_b+1):2*s_b,(s_c+1):2*s_c) = unit_cell;
%  support_abc((s_a+1):2*s_a,(s_b+1):2*s_b,(s_c+1):2*s_c) = ones(s_a,s_b,s_c);
%  save('obj_triple_abc.mat','obj_triple_abc');
%  save('support_abc.mat','support_abc');
%  
%  
%  % Triple in c direction
 obj_triple_c = zeros(s_a,s_b,3*s_c);
 support_c = zeros(s_a,s_b,3*s_c);
 obj_triple_c(:,:,(s_c+1):2*s_c) = unit_cell_abs;
 support_c(:,:,(s_c+1):2*s_c) = ones(s_a,s_b,s_c);
 save('obj_triple_3A.mat','obj_triple_c');
 save('support_triple_3A.mat','support_c');
 
%  % Triple in a direction
%  obj_triple_a = zeros(s_a,s_b,3*s_c);
%  support_a = zeros(s_a,s_b,3*s_c);
%  obj_triple_a((s_a+1):2*s_a,:,:) = unit_cell;
%  support_a((s_a+1):2*s_a,:,:) = ones(s_a,s_b,s_c);
%  save('support_a.mat','support_a');
%  
%  % Triple in b direction
%  obj_triple_b = zeros(s_a,s_b,3*s_c);
%  support_b = zeros(s_a,s_b,3*s_c);
%  obj_triple_b(:,(s_b+1):2*s_b,:) = unit_cell;
%  support_b(:,(s_b+1):2*s_b,:) = ones(s_a,s_b,s_c);
%  save('obj_triple_b.mat','obj_triple_b');
%  save('support_b.mat','support_b');
