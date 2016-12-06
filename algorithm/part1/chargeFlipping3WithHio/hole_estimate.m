% This piece of code is used to find hole in a-b plane.
% Then a support will be imposed based on the hole in center.
% First created by Yun Zhao on Nov 24, 2014
% Last editted by Yun Zhao on Nov 24, 2014
clear all

prompts = 'please give the unit cell file: \n';
unit_cell_file_name = input(prompts,'s');
if isempty(unit_cell_file_name)
    unit_cell_file_name= 'unit_cell_shift_3A.mat';
end
rho_xyz = importdata(unit_cell_file_name);
[s_a,s_b,s_c] = size(rho_xyz);

% Now make a triple cell
triple_cell_shift = zeros(s_a,s_b,3*s_c-2);
triple_cell_shift(:,:,s_c:2*s_c-1) = rho_xyz;
save('triple_cell_shift_3A.mat','triple_cell_shift');



rho_xy = sum(rho_xyz,3);
% imagesc(abs(rho_xy))

% Now make a support.
prompts_2 = 'how many pixels do you want to set as zeros? \n';
solvent_fraction = input(prompts_2);
rho_xy_1d = reshape(rho_xy,[1 s_a*s_b]);
rho_xy_1d_sort = sort(rho_xy_1d);
rho_cutoff = rho_xy_1d_sort(int32(solvent_fraction*s_a*s_b));
support_hole_2d = zeros(s_a,s_b);


for i = 1:s_a
    for j = 1:s_b
        if rho_xy(i,j)>rho_cutoff
            support_hole_2d(i,j) = 1;
        else
            support_hole_2d(i,j) = 0;
        end
    end
end

% Now duplicate it into triple cell.
support_hole_3d = repmat(support_hole_2d,1,1,3*s_c-2);

% Compact support from z
support_z = zeros(s_a,s_b,3*s_c-2);
support_z(:,:,s_c:2*s_c-1)=ones(s_a,s_b,s_c);

% Support for HIO algorithm. Which is a combination of compact support
% along z and hole in the center.
support_hio = support_hole_3d.*support_z;
save('triple_support_hio_3A.mat','support_hio');

% support_file_name = 'support_hio_w_'+int2str(int16(solvent_fraction*100))+'.mat';
% save(support_file_name,'support_hio');

% support_large = zeros(s_a+5,s_b+5,3*s_c+3);
% 
% for i = 1:s_a
%     for j = 1:s_b
%         for k = 1:3*s_c-2
%             support_large(2+i,2+j,2+k)=support_hio(i,j,k);
%         end
%     end
% end
% % support_large(2:(s_a+2),2:(s_b+2),2:(3*s_c)) = support_hio;
% save('support_hio_large.mat','support_large');



