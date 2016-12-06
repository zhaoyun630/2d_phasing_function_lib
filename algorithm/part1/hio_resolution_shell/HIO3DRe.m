%This implementation of the Fienup HIO algorithm in 3D works for real objects. JCHS.


%
% UW/JCHS/YZ.  2014.  Works in Matlab R2014
% Last editted by Yun Zhao on April 7th, 2014

% 3D HiO for real object (feb 03). 
% Read J.Fienup  Applied Optics 21, p.2758 (1982) and  J. Opt Soc Am. A4  p.118 (1987) to understand this program.
% Fienup hio algorithm, first written by Uwe Weierstall. TWO DIMENSIONAL.
% Save experimental modulus in DM as text file. Load it in Matlab and call it "modulus".
% e.g.: modulus = load('DMfilename'); 
% Similarly load support image and call it "support". This program has the modulus constraint in 
% real space. So it assumes a real function for the "object". Hence it would work for a weak phase object.
% 
% obj = load('DM S*Siobj');
%obj=abs(obj);
%save obj
obj=importdata('obj_lysozym_sigma_pdb.mat');   % This is a 512 x 512 two-dimensional array forming the test object. For electron diffraction it is the projected potential
% of a weak phase object. For xrays it is the projected charge density.
modulus1=abs(fftn(obj));     %Discard phases of transform of test object. We will try to find these using the Fienup HIO algorithm.
%modulus = poissrnd(modulus1);
modulus = modulus1;
%support=load('Suport');
%save support
% modulus = importdata('difvol_3.mat');
support = importdata('support_41_41_138.mat');   % This is a mask, which sets to zero intensity outside the approximately known boundary of the object.
phase_prior = importdata('phase_prior.mat');
phase_grid = importdata('phase_grid.mat');
%imshow(obj,[])  % display object
%pause
%Difpat=fftshift(modulus);  % move origin
%imshow(Difpat,[])  % display diffraction pattern intensity.
%pause
%imshow(support,[])  % display support.
%pause
phase_true=angle(fftn(obj));
phi = rand(41,41,138)*2*pi;   % first estimate of phases is random numbers.
% phi = phi.*phase_grid + phase_prior;
%phi = phase_true;
%phi = poissrnd(phase_true);
%phi = zeros(53,53,73);
phase = exp(1i*phi);
phase_error_denominator=sum(sum(sum(modulus.^2)));
%phase = importdata('phase_3.mat');
Gkuv = modulus.*phase;  % initial diff pattern - known amps and random phases.
gprime_kxy = ifftn(Gkuv); % first estimate of image/
gkxy = gprime_kxy.* support;  % apply known support
% imshow(abs(gkxy),[]);
support2 = (support-1)*(-1);  % complement of support.

Gkuv = fftn(gkxy);
Gprime_kuv = modulus.*exp(1i*angle(Gkuv));  % replace with known fourier modulus (diffracted intensities)
% Gprime_kuv = modulus.*exp(1i*(angle(Gkuv).*phase_grid +phase_prior));
gprime_kxy = ifftn(Gprime_kuv);    % next extimate of image.
gk = (abs(gprime_kxy)).^2;  		% error calcn.  Note: DM program does not have square here. now same as 3D
rms1 = sum(sum(sum(gk)));				%first iteration to get rms1 (root mean square of support area after first iteration). 
gkxy = gprime_kxy.*support;			%rms1 is used in subsequent iterations as normalization factor. 
% imshow(abs(gkxy),[]);
   phase_error_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
   phase_error_perc = phase_error_numerator/phase_error_denominator;
   error_info = [1 phase_error_perc];
% fprintf('HiO rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
for w=1:10
 for j = 1:20                      % Do HiO.   change to desired value
   Gkuv = fftn(gkxy);
   Gprime_kuv = modulus.*exp(1i*angle(Gkuv));
   % Gprime_kuv = modulus.*exp(1i*(angle(Gkuv).*phase_grid +phase_prior));
   gprime_kxy = ifftn(Gprime_kuv);
   gk = (abs(gprime_kxy.*support2)).^2;    % error calcn.  DM program does not have square.
   rms = sqrt(sum(sum(sum(gk)))/rms1);  % i updated this line by adding sqrt
   phase_error_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
   phase_error_perc = phase_error_numerator/phase_error_denominator;
   error_info_iter = [rms phase_error_perc];
   error_info = [error_info; error_info_iter];
   fprintf('HiO rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
   gkoutnew = gkxy.* support2 -0.8 * gprime_kxy.* support2;   % feedback factor 0.8
   gkxy = abs(gprime_kxy).*support;     % modulus constraint in real space
   gkxy = gkxy + gkoutnew;
   %b1_1=abs(gkxy);          % added by Yun
   %b1=['w_' int2str(w) 'j_' int2str(j)];   % added by Yun
   %save(b1,'b1_1');         % added by Yun
   %imshow(abs(gkxy),[]);
   %pause(0.1);
 end

  for t = 1:20  % Do Error reduction 10 times.
      Gkuv = fftn(gkxy);
      Gprime_kuv = modulus.*exp(1i*angle(Gkuv));  % replace with known fourier modulus (diffracted intensities)
      % Gprime_kuv = modulus.*exp(1i*(angle(Gkuv).*phase_grid +phase_prior));
   	gprime_kxy = ifftn(Gprime_kuv);
   	gk = (abs(gprime_kxy.*support2)).^2;  
   	rms = sqrt(sum(sum(sum(gk)))/rms1);
    phase_error_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
    phase_error_perc = phase_error_numerator/phase_error_denominator;
    error_info_iter = [rms phase_error_perc];
    error_info = [error_info; error_info_iter];
    fprintf('error-reduction rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
   	%fprintf('error-reduction: rms = %g.\n',rms);
      gkxy = abs(gprime_kxy).*support;
	%imshow(abs(gkxy),[]);
   end
end
%imshow(abs(gkxy),[]);
save('error_info.mat','error_info');
