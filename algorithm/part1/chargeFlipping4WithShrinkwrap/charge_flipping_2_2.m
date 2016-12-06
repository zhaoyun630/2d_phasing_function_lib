% Ab initio phasing with charge flipping algorithm
% This piece of code exclusively use charge flipping
% First created by Yun Zhao on Nov 12th, 2014
% Last editted by Yun Zhao on Jan 18th, 2015
% This is the edition I wrote at Xiantao, China

% Compared with version 2_1, this version takes hkl file as input.


% obj=importdata('unit_cell_shift.mat');   
% F_modulus_shifted_observed=abs(fftn(obj));     %Discard phases of transform of test object. We will try to find these using the Fienup HIO algorithm.
% % after fftshift, the center is F(000).
% %  F_modulus_centered_observed = fftshift(modulus);

F_mask_centered = importdata('rec_lattice_mask_1A.mat');
F_mask_uncentered = ifftshift(F_mask_centered);

F_complex_centered_observed = importdata('rec_lattice_1A.mat');
F_modulus_centered_observed = abs(F_complex_centered_observed);
F_modulus_shifted_observed = ifftshift(F_modulus_centered_observed);
obj = real(ifft(ifftshift(F_complex_centered_observed)));


% set some global variables.
solvent_level = 0.7;
flip_frac = 0.79;
weak_amp_percent = 0.001;
iter_num_total = 5000;
R_factor = [];

[s_a,s_b,s_c] = size(obj);
flip_index = floor(flip_frac*s_a*s_b*s_c);

% Check if s_a,s_b,s_c are even or odd.
% Actually, they all should be odd. As h,k,l will be taken with values:
% 0,+/-1,+/-2,....
s_a_group = mod(s_a,2);
s_b_group = mod(s_b,2);
s_c_group = mod(s_c,2);

% The matrix index at the center.
s_a_center = (s_a+s_a_group)/2;
s_b_center = (s_b+s_b_group)/2;
s_c_center = (s_c+s_c_group)/2;

% Creat random phase with centrosymmetry.
phase_0 = 2*pi*rand(s_a,s_b,s_c);

% The phase of (h00),(0k0),(00l) should be zero.
% phase_0(s_a_center,:,:)=0;
% phase_0(:,s_b_center,:)=0;
% phase_0(:,:,s_c_center)=0;
phase_0(s_a_center,s_b_center,s_c_center)=0;

% Set our starting value F(000) = 0.
% In later iterations, it will be kept as its value in previous iteration.
modulus_centered_0 = F_modulus_centered_observed;
modulus_centered_0(s_a_center,s_b_center,s_c_center)=0;

% modulus_centered_0 = fftshift(F_modulus_shifted_observed);
% modulus_centered_0(s_a_center,s_b_center,s_c_center)=0;

% Now apply Friedel's symmetry to our starting phase:phi(h,k,l) = - phi(h,k,l);
for i = 1: (s_a+s_a_group)/2
    for j = 1:s_b
        for k = 1:s_c
            phase_0(i,j,k) = - phase_0(s_a-i+1,s_b-j+1,s_c-k+1);
        end
    end
end

% First estimation of our object.
F_comp_0 = modulus_centered_0.*exp(1i*phase_0);
% obj_comp = ifftn(F_comp_0);
% obj_k = abs(obj_comp); % First estimate of object
obj_k = real(ifftn(ifftshift(F_comp_0)));

% Start Charge flipping iteration
obj_1d_array = reshape(obj_k,[1,s_a*s_b*s_c]);
rho_sigma = std(obj_1d_array);
rho_mean = mean(obj_1d_array);
rho_max = max(obj_1d_array);


% rho_cut_off = solvent_level*rho_sigma;
rho_sort = sort(obj_1d_array);
% rho_cut_off = rho_sort(flip_index);
% rho_cut_off = 0.01;

for iter = 1: iter_num_total
    % Charge flipping in real space
    obj_1d_array = reshape(obj_k,[1,s_a*s_b*s_c]);
    rho_sigma = std(obj_1d_array);
    rho_mean = mean(obj_1d_array);
    rho_max = max(obj_1d_array);
    rho_sort = sort(obj_1d_array);
    rho_cut_off = rho_sort(flip_index);
%     rho_cut_off = solvent_level*rho_sigma;
    obj_k_1= obj_k;
    for i = 1:s_a
        for j = 1:s_b
            for k = 1:s_c
                if abs(obj_k(i,j,k))<rho_cut_off
                    obj_k_1(i,j,k) = -obj_k(i,j,k);
                end
            end
        end
    end
    
    % Constraints in Fourier Space
    F_comp_1 =fftn(obj_k_1);
    F_modulus_1= abs(F_comp_1);
    F_modulus_constraint = F_modulus_shifted_observed;
    F_modulus = F_modulus_1.*F_mask_uncentered;
    r_i = sum(sum(sum(abs(abs(F_modulus)-abs(F_modulus_constraint)))))/(sum(sum(sum(abs(F_modulus_constraint))))-abs(F_modulus_constraint(1,1,1)));
    R_factor = [R_factor;r_i];
    phase = angle(F_comp_1);
    F_modulus_constraint(1,1,1)=F_modulus_1(1,1,1);
    
%     % Find the structure amplitudes below than certain percentage.
%     F_modulus_1d = reshape(F_modulus,[1,s_a*s_b*s_c]);
%     F_modulus_1d_sort = sort(F_modulus_1d);
%     weak_index = floor(weak_amp_percent*s_a*s_b*s_c);
%     F_modulus_weak = F_modulus_1d_sort(weak_index);
%     
%     % Replace phase associated with weak Bragg spots to pi/2
%     for i = 1:s_a
%         for j = 1:s_b
%             for k = 1:s_c
%                 if F_modulus(i,j,k)<F_modulus_weak
%                     phase(i,j,k) = phase(i,j,k)+pi/2;
%                 end
%             end
%         end
%     end
    
    
    % Keep the F(000) from the last iteration
    % Note that F(000) does not locate at the center of this 3D array.
    % It is the first in FFT from Matlab
    F_comp_2 = F_modulus_constraint.*exp(1i*phase);
    obj_k = real(ifftn(F_comp_2));
    
end

save('obj_k_5000iter_hkl_1A_sf79.mat','obj_k');
save('r_factor_5000iter_hkl_1A_sf79.mat','R_factor');

R_factor = importdata('r_factor_5000iter_hkl_1A_sf79.mat');
iter_num_total = length(R_factor);
iter_num = linspace(1,iter_num_total,iter_num_total);
plt1 = figure(1);
plot(iter_num,R_factor,'-b.');
lg1_text1 = strcat('total iteration: ',int2str(iter_num_total)); 
lg1_text2 = strcat('weak reflection percentage:' , int2str(weak_amp_percent));
% legend(lg1_text1,lg1_text2);
% text(700,1.54,lg1_text1,'FontSize', 15,'FontWeight','bold');
% text(700,1.53,lg1_text2,'FontSize', 15,'FontWeight','bold');
xlabel('iteration','FontSize', 15,'FontWeight','bold');
ylabel('R factor','FontSize', 15,'FontWeight','bold');
% rectangle('Position',[690,1.525,250,0.022]);
saveas(plt1,'r_factor_hkl_iteration_5000_1A_sf79.tif');






