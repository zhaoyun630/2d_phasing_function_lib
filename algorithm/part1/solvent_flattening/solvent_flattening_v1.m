% Solvent flattening with a molecular envelope constraint.

%
% First Created by Yun Zhao on April 17th, 2014
% Last editted by Yun Zhao on Feb 10th, 2015
% Compare with Version 2, it takes real value off inverse FFT from structure factors
clear all

obj_model=abs(importdata('obj_abs_triple_3A.mat'));   
% obj_model=abs(importdata('obj_triple_c_sigma_2.mat'));  
F_modulus_shifted_observed_1 = abs(fftn(obj_model));     %Discard phases of transform of test object. We will try to find these using the Fienup HIO algorithm.
%modulus = poissrnd(modulus1);
F_modulus_shifted_observed = F_modulus_shifted_observed_1;
support_1 = importdata('support_shape_pc_25.mat');   % This is a mask, which sets to zero intensity outside the approximately known boundary of the object.
[s_a,s_b,s_c] = size(support_1);
rms_cc_Rfactor = [];

F_angle_true=angle(fftn(obj_model));
F_angle_initial = rand(s_a,s_b,s_c)*2*pi;   % first estimate of phases is random numbers.




F_phase_exp_term = exp(1j*F_angle_initial);
phase_cc_denominator=sum(sum(sum(F_modulus_shifted_observed.^2)));
%phase = importdata('phase_3.mat');
F_complex_shifted_prime = F_modulus_shifted_observed.*F_phase_exp_term;  % initial diff pattern - known amps and random phases.
obj_estimate_prime = ifftn(F_complex_shifted_prime); % first estimate of image/
obj_estimate = obj_estimate_prime.* support_1;  % apply known support
% imshow(abs(gkxy),[]);
support_2 = 1-support_1;  % complement of support.

F_complex_shifted_prime = fftn(obj_estimate);
%Gprime_kuv = modulus.*exp(1i*angle(Gkuv));  % replace with known fourier modulus (diffracted intensities)
 F_complex_shifted_modulus_constraint = F_modulus_shifted_observed.*exp(1i*(angle(F_complex_shifted_prime)));
obj_estimate_prime = ifftn(F_complex_shifted_modulus_constraint);    % next extimate of image.
obj_outside_support = (abs(obj_estimate_prime.*support_2)).^2;  		% error calcn.  Note: DM program does not have square here. now same as 3D
rms = sqrt(sum(sum(sum(obj_outside_support)))/sum(sum(sum(abs(obj_estimate_prime)))));				%first iteration to get rms1 (root mean square of support area after first iteration). 
obj_estimate = obj_estimate_prime.*support_1;			%rms1 is used in subsequent iterations as normalization factor. 

phase_cc_numerator=sum(sum(sum(F_modulus_shifted_observed.^2.*cos(angle(F_complex_shifted_prime)-F_angle_true))));
%   phase_cc_numerator=sum(sum(sum(modulus.^2.*abs(cos(angle(Gkuv)-phase_true)))));
phase_cc = phase_cc_numerator/phase_cc_denominator;



F_modulus = abs(F_complex_shifted_prime);
F_modulus_constraint = F_modulus_shifted_observed;
% F_modulus_constraint(1,1,1)=F_modulus(1,1,1);
F_modulus(1,1,1) = F_modulus_constraint(1,1,1);
r_factor = sum(sum(sum(abs(abs(F_modulus)-abs(F_modulus_constraint)))))/(sum(sum(sum(abs(F_modulus_constraint))))-abs(F_modulus_constraint(1,1,1)));

rms_cc_Rfactor = [rms_cc_Rfactor; rms phase_cc r_factor];
  
  
% fprintf('HiO rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
for w=1:15
    for j = 1:40                      % Do HiO.   change to desired value
       F_complex_shifted_prime = fftn(obj_estimate);
       F_complex_shifted_modulus_constraint = F_modulus_shifted_observed.*exp(1i*angle(F_complex_shifted_prime));
       obj_estimate_prime = ifftn(F_complex_shifted_modulus_constraint);
       obj_outside_support = (abs(obj_estimate_prime.*support_2)).^2;    % error calcn.  DM program does not have square.
       rms = sqrt(sum(sum(sum(obj_outside_support)))/sum(sum(sum(abs(obj_estimate_prime)))));  % i updated this line by adding sqrt
       phase_cc_numerator=sum(sum(sum(F_modulus_shifted_observed.^2.*cos(angle(F_complex_shifted_prime)-F_angle_true))));
%        phase_cc_numerator=sum(sum(sum(modulus.^2.*abs(cos(angle(Gkuv)-phase_true)))));
       phase_cc = phase_cc_numerator/phase_cc_denominator;
       F_modulus = abs(F_complex_shifted_prime);
       F_modulus_constraint = F_modulus_shifted_observed;
%        F_modulus_constraint(1,1,1)=F_modulus(1,1,1);
       F_modulus(1,1,1) = F_modulus_constraint(1,1,1);
       r_factor = sum(sum(sum(abs(abs(F_modulus)-abs(F_modulus_constraint)))))/(sum(sum(sum(abs(F_modulus_constraint))))-abs(F_modulus_constraint(1,1,1)));
       rms_cc_info_iter = [rms phase_cc r_factor];
       rms_cc_Rfactor = [rms_cc_Rfactor; rms_cc_info_iter];
       % fprintf('HiO rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
       gkoutnew = obj_estimate.* support_2 -0.8 * obj_estimate_prime.* support_2;   % feedback factor 0._obj_0sigma_sf_008
       obj_estimate_1 = abs(obj_estimate_prime).*support_1;     % modulus constraint in real space
       obj_estimate = obj_estimate_1 + gkoutnew;
       %pause(0.1);
    end

    for t = 1:20  % Do Error reduction 10 times.
        F_complex_shifted_prime = fftn(obj_estimate);
        % Gprime_kuv = modulus.*exp(1i*angle(Gkuv));  % replace with known fourier modulus (diffracted intensities)
        F_complex_shifted_modulus_constraint = F_modulus_shifted_observed.*exp(1i*(angle(F_complex_shifted_prime)));
    	obj_estimate_prime = ifftn(F_complex_shifted_modulus_constraint);
    	obj_outside_support = (abs(obj_estimate_prime.*support_2)).^2;  
        rms = sqrt(sum(sum(sum(obj_outside_support)))/sum(sum(sum(abs(obj_estimate_prime)))));
        phase_cc_numerator=sum(sum(sum(F_modulus_shifted_observed.^2.*cos(angle(F_complex_shifted_prime)-F_angle_true))));
        phase_cc = phase_cc_numerator/phase_cc_denominator;
        F_modulus = abs(F_complex_shifted_prime);
        F_modulus_constraint = F_modulus_shifted_observed;
%         F_modulus_constraint(1,1,1)=F_modulus(1,1,1);
        F_modulus(1,1,1) = F_modulus_constraint(1,1,1);
        r_factor = sum(sum(sum(abs(abs(F_modulus)-abs(F_modulus_constraint)))))/(sum(sum(sum(abs(F_modulus_constraint))))-abs(F_modulus_constraint(1,1,1)));
        rms_cc_info_iter = [rms phase_cc r_factor];
        rms_cc_Rfactor = [rms_cc_Rfactor; rms_cc_info_iter];
        % fprintf('error-reduction rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
    	% fprintf('error-reduction: rms = %g.\n',rms);
        obj_estimate = abs(obj_estimate_prime).*support_1;
    end
end

% obj = zeros([s_a,s_b,s_c]/3);

% Here support size is triple as the original one in each dimension.
% We just take the middle part, which is reconstructed obj.
obj_hio = obj_estimate(:,:,int32(s_c/3+1):int32(s_c*2/3)); 

save('obj_hio_obj_0sigma_sf_008.mat','obj_hio');
% save('gkxy_sf.mat','gkxy');
save('rms_hio_obj_0sigma_sf_008.mat','rms_cc_Rfactor');

load('rms_hio_obj_0sigma_sf_008.mat');
figure(1);
% rms_cc = importdata('rms_hio_hole_5.mat');
L = size(rms_cc_Rfactor);
x = linspace(1,L(1),L(1));
plt1=plot(x,rms_cc_Rfactor(:,1),'b-*',x,rms_cc_Rfactor(:,2),'r-x');
x1=xlabel('interation');
set(x1,'FontSize', 15);
set(x1,'FontWeight','bold');
y1 = ylabel('CC and rms');
set(y1,'FontSize',15);
set(y1,'FontWeight','bold')
fleg = legend('rms','CC');
set(fleg,'FontSize',15);
set(fleg,'FontWeight','bold');
saveas(gcf,'rms_cc_obj_0sigma_sf_008.tif');


figure(2);
plt2=plot(x,rms_cc_Rfactor(:,3),'b-x');
x1=xlabel('interation');
set(x1,'FontSize', 15);
set(x1,'FontWeight','bold');
y1 = ylabel('R_factor');
set(y1,'FontSize',15);
set(y1,'FontWeight','bold')
% fleg = legend('rms','CC','R factor');
% set(fleg,'FontSize',15);
% set(fleg,'FontWeight','bold');
saveas(plt2,'rfactor_obj_0sigma_sf_008.tif');

