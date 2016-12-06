% Ab initio phasing with charge flipping algorithm
% This piece of code exclusively use charge flipping
% First created by Yun Zhao on Nov 12th, 2014
% Last editted by Yun Zhao on Nov 12th, 2014

clear all

obj=importdata('unit_cell_1A.mat');   
modulus=abs(fftn(obj));     %Discard phases of transform of test object. We will try to find these using the Fienup HIO algorithm.

% set some global variables.
solvent_level = 0.9;
weak_amp_percent = 0.3;
iter_num_total = 400;
R_factor = [];
% after fftshift, the center is F(000).
modulus_centered = fftshift(modulus);

[s_a,s_b,s_c] = size(obj);

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
modulus_centered_0 = modulus_centered;
modulus_centered_0(s_a_center,s_b_center,s_c_center)=0;

% Now apply Friedel's symmetry to our starting phase:phi(h,k,l) = - phi(h,k,l);
for i = 1: (s_a-s_a_group)/2
    for j = 1:s_b
        for k = 1:s_c
            phase_0(i,j,k) = - phase_0(s_a-i+1,s_b-j+1,s_c-k+1);
        end
    end
end

% First estimation of our object.
F_comp_0 = modulus_centered_0.*exp(1i*phase_0);
obj_comp = ifftn(ifftshift(F_comp_0));
obj_k = real(obj_comp); % First estimate of object
obj_k_1 = obj_k;


% Start Charge flipping iteration
obj_1d_array = reshape(obj_k,[1,s_a*s_b*s_c]);
rho_sigma = std(obj_1d_array);
rho_mean = mean(obj_1d_array);
rho_max = max(obj_1d_array);

rho_cut_off = solvent_level*rho_sigma;

for iter = 1: iter_num_total
    % Charge flipping in real space
    obj_1d_array = reshape(obj_k,[1,s_a*s_b*s_c]);
    rho_sigma = std(obj_1d_array);
    rho_mean = mean(obj_1d_array);
    rho_max = max(obj_1d_array);
    
    rho_cut_off = solvent_level*rho_sigma;
    for i = 1:s_a
        for j = 1:s_b
            for k = 1:s_c
                if obj_k(i,j,k)<rho_cut_off
                    obj_k_1(i,j,k) = -obj_k(i,j,k);
                end
            end
        end
    end
    
    % Constraints in Fourier Space
    F_comp_1 =fftn(obj_k_1);
    F_modulus = real(F_comp_1);
    phase = angle(F_comp_1);
    R_iter = sum(sum(sum(abs(F_modulus-modulus))))/(sum(sum(sum(abs(modulus))))-modulus(1,1,1));
    R_factor = [R_factor; R_iter];
    
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
%                     phase(i,j,k) = pi/2;
%                 end
%             end
%         end
%     end
    
    
    % Keep the F(000) from the last iteration
    % Note that F(000) does not locate at the center of this 3D array.
    % It is the first in FFT from Matlab
    modulus(1,1,1)=F_comp_1(1,1,1);
    F_comp_2 = modulus.*exp(1i*phase);
    obj_k = real(ifftn(F_comp_2));
    
end

save('obj_k.mat','obj_k');

figure(1)
iter_num = linspace(1,iter_num_total,iter_num_total);
plt1 = plot(iter_num,R_factor);








