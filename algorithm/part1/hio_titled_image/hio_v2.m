%This implementation of the Fienup HIO algorithm in 3D works for real objects. JCHS.


%
% UW/JCHS/YZ.  2014.  Works in Matlab R2014
% Last editted by Yun Zhao on April 17th, 2014
clear all

obj=importdata('obj_lysozym_sigma_pdb.mat');   
modulus1=abs(fftn(obj));     %Discard phases of transform of test object. We will try to find these using the Fienup HIO algorithm.
%modulus = poissrnd(modulus1);
modulus = modulus1;
support = importdata('support_41_41_138.mat');   % This is a mask, which sets to zero intensity outside the approximately known boundary of the object.
phase_prior = importdata('phase_prior.mat');
phase_grid = importdata('phase_grid.mat');

phase_true=angle(fftn(obj));
phi = rand(41,41,138)*2*pi;   % first estimate of phases is random numbers.
 phi = phi.*phase_grid + phase_prior;
% phi = phase_true;
%phi = poissrnd(phase_true);
%phi = zeros(53,53,73);
phase = exp(1i*phi);
phase_cc_denominator=sum(sum(sum(modulus.^2)));
%phase = importdata('phase_3.mat');
Gkuv = modulus.*phase;  % initial diff pattern - known amps and random phases.
gprime_kxy = ifftn(Gkuv); % first estimate of image/
gkxy = gprime_kxy.* support;  % apply known support
% imshow(abs(gkxy),[]);
support2 = (support-1)*(-1);  % complement of support.

Gkuv = fftn(gkxy);
%Gprime_kuv = modulus.*exp(1i*angle(Gkuv));  % replace with known fourier modulus (diffracted intensities)
 Gprime_kuv = modulus.*exp(1i*(angle(Gkuv).*phase_grid +phase_prior));
gprime_kxy = ifftn(Gprime_kuv);    % next extimate of image.
gk = (abs(gprime_kxy)).^2;  		% error calcn.  Note: DM program does not have square here. now same as 3D
rms1 = sum(sum(sum(gk)));				%first iteration to get rms1 (root mean square of support area after first iteration). 
gkxy = gprime_kxy.*support;			%rms1 is used in subsequent iterations as normalization factor. 
% imshow(abs(gkxy),[]);
 %   phase_cc_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
      phase_cc_numerator=sum(sum(sum(modulus.^2.*abs(cos(angle(Gkuv)-phase_true)))));
   phase_cc = phase_cc_numerator/phase_cc_denominator;
   rms_cc = [1 phase_cc];
% fprintf('HiO rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
for w=1:20
 for j = 1:20                      % Do HiO.   change to desired value
   Gkuv = fftn(gkxy);
  % Gprime_kuv = modulus.*exp(1i*angle(Gkuv));
    Gprime_kuv = modulus.*exp(1i*(angle(Gkuv).*phase_grid +phase_prior));
   gprime_kxy = ifftn(Gprime_kuv);
   gk = (abs(gprime_kxy.*support2)).^2;    % error calcn.  DM program does not have square.
   rms = sqrt(sum(sum(sum(gk)))/rms1);  % i updated this line by adding sqrt
   %phase_cc_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
      phase_cc_numerator=sum(sum(sum(modulus.^2.*abs(cos(angle(Gkuv)-phase_true)))));
   phase_cc = phase_cc_numerator/phase_cc_denominator;
   rms_cc_info_iter = [rms phase_cc];
   rms_cc = [rms_cc; rms_cc_info_iter];
  % fprintf('HiO rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
   gkoutnew = gkxy.* support2 -0.8 * gprime_kxy.* support2;   % feedback factor 0.8
   gkxy = abs(gprime_kxy).*support;     % modulus constraint in real space
   gkxy = gkxy + gkoutnew;
   %pause(0.1);
 end

  for t = 1:20  % Do Error reduction 10 times.
      Gkuv = fftn(gkxy);
     % Gprime_kuv = modulus.*exp(1i*angle(Gkuv));  % replace with known fourier modulus (diffracted intensities)
       Gprime_kuv = modulus.*exp(1i*(angle(Gkuv).*phase_grid +phase_prior));
   	gprime_kxy = ifftn(Gprime_kuv);
   	gk = (abs(gprime_kxy.*support2)).^2;  
   	rms = sqrt(sum(sum(sum(gk)))/rms1);
    phase_cc_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
    phase_cc = phase_cc_numerator/phase_cc_denominator;
    rms_cc_info_iter = [rms phase_cc];
    rms_cc = [rms_cc; rms_cc_info_iter];
   % fprintf('error-reduction rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
   	%fprintf('error-reduction: rms = %g.\n',rms);
      gkxy = abs(gprime_kxy).*support;
   end
end

save('gkxy.mat','gkxy');
save('rms_cc.mat','rms_cc');
