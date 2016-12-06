% This piece of code is used to find hole in a-b plane.
% Then a support will be imposed based on the hole in center.
% First created by Yun Zhao on Nov 24, 2014
% Last editted by Yun Zhao on Feb 13, 2014
% This version V3 calculate 3D fraction first, then
% estimate molecular envelope.
% In V1, the 2d projection is calculated first, then estimate hole.

clear all

prompts = 'please give the unit cell file: \n';
unit_cell_file_name = input(prompts,'s');
if isempty(unit_cell_file_name)
    unit_cell_file_name= 'unit_cell_abs_3A.mat';
end
rho_xyz = importdata(unit_cell_file_name);
[s_a,s_b,s_c] = size(rho_xyz);


prompts_2 = 'how many pixels do you want to set as zeros? \n';
solvent_fraction_int = input(prompts_2);
solvent_fraction = solvent_fraction_int/100;
rho_xyz_1d = reshape(rho_xyz,[1 s_a*s_b*s_c]);
rho_xyz_1d_sort = sort(rho_xyz_1d);
rho_cutoff = rho_xyz_1d_sort(int64(solvent_fraction*s_a*s_b*s_c));


% Now obtain a new object after setting certain fraction to zero
rho_solvent_flattening = zeros(s_a,s_b,s_c);
for i = 1:s_a
    for j = 1:s_b
        for k = 1:s_c
            if rho_xyz(i,j,k)>rho_cutoff
                rho_solvent_flattening(i,j,k) = 1;
            else
                rho_solvent_flattening(i,j,k) = 0;
            end
        end
    end
end

% Compact support from z
support_z = zeros(s_a,s_b,3*s_c-2);
support_z(:,:,s_c:2*s_c-1)=rho_solvent_flattening;

% Support for HIO algorithm. Which is a combination of compact support
% along z and hole in the center.
support_hio = support_z;
% save('triple_support_hio_3A_sf_60.mat','support_hio');

support_file_name = strcat('triple_support_hio_3A_sf_',int2str(int32(solvent_fraction_int)),'.mat');
save(support_file_name,'support_hio');

% To make a better 3D view in my show_3d_obj.m 
support_large = zeros(s_a+3,s_b+3,3*s_c+1);

for i = 1:s_a
    for j = 1:s_b
        for k = 1:3*s_c-2
            support_large(2+i,2+j,2+k)=support_hio(i,j,k);
        end
    end
end
% support_large(2:(s_a+2),2:(s_b+2),2:(3*s_c)) = support_hio;
save('support_hio_large.mat','support_large');


