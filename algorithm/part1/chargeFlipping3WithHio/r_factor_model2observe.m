% This piece of code is used to calculate R_factor between model
% and experimental structure factors.

 prompts = ' Please give a model file: \n';
 model_file_name = input(prompts,'s');
 if isempty(model_file_name)
     model_file_name = 'unit_cell_abs.mat';
 end
 
 obj = importdata(model_file_name);
 
 prompts2 = ' Please give a model file: \n';
 rho_shift_file_name = input(prompts2,'s');
 if isempty(rho_shift_file_name)
%      rho_shift = importdata('rho_shift.mat');
     rho_shift = 0;
 end
 
 F_mask_centered = importdata('rec_lattice_mask.mat');
 F_mask_uncentered = ifftshift(F_mask_centered);

%  rho_shift = 0;
 [s_a,s_b,s_c]=size(obj);
 const_shift_obj = rho_shift*ones(s_a,s_b,s_c);
 
 % Complex structure factor from given structure
 F_obj_comp_shift = fft(obj); 
%  F_const_shift_obj = fft(const_shift_obj);
 F_const_shift_obj = 0;
 F_obj_comp_original = F_obj_comp_shift + F_const_shift_obj;
 
 % Complex structure factor (observed amplitude + phase from model)
 F_obs_comp = importdata('rec_lattice.mat');
 
 F_obj_amplitude = abs(F_obj_comp_original);
 F_obs_amplitude = ifftshift(abs(F_obs_comp));
 
 % Mask the F_000
 F_obj_amplitude_without_center = F_obj_amplitude;
 F_obj_amplitude_without_center(1,1,1) = 0;
 
 F_obs_amplitude_without_center = F_obs_amplitude;
 F_obs_amplitude_without_center(1,1,1) = 0;
 
 save('F_obj_amplitude_without_center.mat','F_obj_amplitude_without_center');
 save('F_obs_amplitude_without_center.mat','F_obs_amplitude_without_center');
 
 r_i = sum(sum(sum(abs(abs(F_obj_amplitude_without_center)-abs(F_obs_amplitude_without_center)))))/sum(sum(sum(abs(F_obs_amplitude_without_center))));
 
 
 fprintf('R factor is: %f \n', r_i);
 
 
 