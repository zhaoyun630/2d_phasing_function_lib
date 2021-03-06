% initial unit cell

% F = difvol_1.*phase_1;
% obj_complex = ifftn(F);
% obj_complex_centered = fftshift(obj_complex);
% obj_c = abs(obj_complex_centered);
% obj_unc = abs(obj_complex);
% obj_c1 = real(obj_complex_centered);
% obj_c2 = imag(obj_complex_centered);
% 
% % Tripled unit cell
% 
% F_tri = difvol_3abc.*phase_3abc;
% tri_obj_complex = ifftn(F_tri);
% tri_obj_complex_centered = fftshift(tri_obj_complex);
% tri_obj_c = abs(tri_obj_complex_centered);
% tri_obj_unc = abs(tri_obj_complex);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial unit cell

F = difvol_1_full.*phase_1_full;
obj_complex = ifftn(F);
obj_complex_centered = fftshift(obj_complex);
obj_c = abs(obj_complex_centered);
obj_unc = abs(obj_complex);
obj_c1 = real(obj_complex_centered);
obj_c2 = imag(obj_complex_centered);

% Tripled unit cell

F_tri = difvol_3abc_full.*phase_3abc_full;
tri_obj_complex = ifftn(F_tri);
tri_obj_complex_centered = fftshift(tri_obj_complex);
tri_obj_c = abs(tri_obj_complex_centered);
tri_obj_unc = abs(tri_obj_complex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now insert two zero planes along z direction
% difvol_nc = zeros(a,b,3*c-2);
% phase_nc = ones(a,b,3*c-2);
% for i = 1:c
%     num = 3*i-2; % the ith row in Difvol corresponds to 1+3(i-1) in Difvol
%     difvol_nc(:,:,num)=difvol_1(:,:,i);
%     phase_nc(:,:,num) = phase_1(:,:,i);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show difference between obj and obj_complex in fourier transform
% diff_vol = fftn(obj_unc);
% diff_intensity1 = abs(diff_vol);
% diff_angle1 = angle(diff_vol);
% diff_phase = exp(1i*diff_angle);
% 
% abs_diff_phase = abs(exp(1i*(diff_angle1(20,20,18:38))))-abs(exp(1i*angle_1(20,20,18:38)));
% 
% abs_diff_intensity_perc = abs(diff_intensity1(20,20,18:38) - difvol_1(20,20,18:38))./(diff_intensity1(20,20,18:38) + difvol_1(20,20,18:38));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the standard deviation
% gkxy_in = importdata('gkxy_9.mat');
gkxy_in = importdata('unit_cell.mat');
gkxy_ave = mean(mean(mean(gkxy_in)));
gkxy_std = sqrt( mean(mean(mean((gkxy_in- gkxy_ave).^2))));





