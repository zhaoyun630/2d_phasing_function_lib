% This piece of code is used to calcuate the R_factor between model and
% the model after solvent flattening approaximation.

clear all
prompts = 'Please give the obj file: \n';
obj_file = input(prompts,'s');
if isempty(obj_file)
    obj_file = 'unit_cell_1A.mat'; % default file, if not input is given
%     obj_file = 'obj_k.mat';
end
obj = importdata(obj_file);
[s_a,s_b,s_c] = size(obj);

prompts = 'Please set the sigma level: \n';
sigma_level = input(prompts);
if isempty(sigma_level)
    sigma_level = 3.3;
end

rho_1D_array = reshape(obj,1,s_a*s_b*s_c); % Squeeze 3d array to 1d array
rho_mean = mean(rho_1D_array);
rho_sigma = std(rho_1D_array);
rho_cut_off = sigma_level*rho_sigma;
rho_cut_off = 0.05;


obj_sovent_flattening= obj;
obj_charge_flipping = obj;
for i = 1:s_a
    for j = 1:s_b
        for k = 1:s_c
            if obj(i,j,k)<rho_cut_off
                obj_sovent_flattening(i,j,k) = 0;
                obj_charge_flipping(i,j,k) = -obj(i,j,k);
            end
        end
    end
end

F_obj = abs(fft(obj));
F_obj_sf = abs(fft(obj_sovent_flattening));
F_obj_cf = abs(fft(obj_charge_flipping));

% Calculate R-factor
R_factor_sf = sum(sum(sum(abs(F_obj_sf-F_obj))))/sum(sum(sum(abs(F_obj))));
R_factor_cf = sum(sum(sum(abs(F_obj_cf-F_obj))))/sum(sum(sum(abs(F_obj))));

% Calculate R-factor value without F(000)
F_obj(1,1,1)=0;
F_obj_sf(1,1,1)= 0;
F_obj_cf(1,1,1)= 0;
R_factor_sf_without_center = sum(sum(sum(abs(F_obj_sf-F_obj))))/sum(sum(sum(abs(F_obj))));
R_factor_cf_without_center = sum(sum(sum(abs(F_obj_cf-F_obj))))/sum(sum(sum(abs(F_obj))));

fprintf('For the model with solvent flatting, the R_factor is %f \n',R_factor_sf);
fprintf(' the R_factor without considering F(000) is %f \n',R_factor_sf_without_center);
fprintf('For the model with charge flipping, the R_factor is %f \n',R_factor_cf);
fprintf(' the R_factor without considering F(000) is %f \n',R_factor_cf_without_center);

sigma_level_int = int32(sigma_level);
filename1 = strcat('obj_sf_sigma_',int2str(sigma_level_int),'.mat');
filename2 = strcat('obj_cf_sigma_',int2str(sigma_level_int),'.mat');
save(filename1,'obj_sovent_flattening');
save(filename2,'obj_charge_flipping');







