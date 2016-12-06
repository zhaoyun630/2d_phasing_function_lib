% Estimate a support from obj.
% The density will be zero if p < p_cutoff.


fprintf(' make sure the input obj file is a triple cell along c \n');
fprintf(' otherwise, it wouldn''t work properly. \n')
prompts = 'please give the obj file: \n';
obj_file = input(prompts,'s');
if isempty(obj_file)
    obj_file = 'obj_triple_3A.mat';
end
obj = importdata(obj_file);
[s_a,s_b,s_c] = size(obj);

notes1 = 'please set the sigma level for obj: \n';
sigma_level_obj = input(notes1);

notes2 = 'please set the sigma level for support: \n';
sigma_level_support = input(notes2);

obj_middle = obj(:,:,s_c/3+1:s_c*2/3);
rho_1D_array = reshape(obj_middle,1,s_a*s_b*s_c/3); % Squeeze 3d array to 1d array
rho_mean = mean(rho_1D_array);
rho_sigma = std(rho_1D_array);

support = zeros(s_a,s_b,s_c);
obj_sigma = zeros(s_a,s_b,s_c);

rho_cutoff_obj = sigma_level_obj*rho_sigma;
rho_cutoff_support = sigma_level_support*rho_sigma;

% Get the new object.
for i=1:s_a
    for j = 1:s_b
        for k=1:s_c
            if obj(i,j,k) < rho_cutoff_obj
                obj_sigma(i,j,k) = 0;
            else
                obj_sigma(i,j,k) = obj(i,j,k);
            end
        end
    end
end

% Get the new support
for i=1:s_a
    for j = 1:s_b
        for k=1:s_c
            if obj(i,j,k) < rho_cutoff_support
                support(i,j,k) = 0;                
            else
                support(i,j,k) = 1;
            end
        end
    end
end


% obj_sigma_file = 'obj_triple_c_sigma_' + int2str(sigma_level)+'.mat';
obj_sigma_file = strcat('obj_triple_c_sigma_2','.mat');
save('support_sigma_0_7.mat', 'support');
save(obj_sigma_file,'obj_sigma');
