%This implementation of the Fienup HIO algorithm in 2D works for real objects. JCHS.


%
% UW/JCHS.  2003.  June 2010. Works in Matlab R2009

%2D HiO for real object (feb 03). 
% Read J.Fienup  Applied Optics 21, p.2758 (1982) and  J. Opt Soc Am. A4  p.118 (1987) to understand this program.
% Fienup hio algorithm, first written by Uwe Weierstall. TWO DIMENSIONAL.
% Save experimental modulus in DM as text file. Load it in Matlab and call it "modulus".
% e.g.: modulus = load('DMfilename'); 
% Similarly load support image and call it "support". This program has the modulus constraint in 
% real space. So it assumes a real function for the "object". Hence it would work for a weak phase object.
% 
%obj = load('DM S*Siobj');
%obj=abs(obj);
%save obj
load obj.mat   % This is a 512 x 512 two-dimensional array forming the test object. For electron diffraction it is the projected potential
% of a weak phase object. For xrays it is the projected charge density.
modulus=abs(fft2(obj));     %Discard phases of transform of test object. We will try to find these using the Fienup HIO algorithm.
%support=load('Suport');
%save support
support = importdata('support');   % This is a mask, which sets to zero intensity outside the approximately known boundary of the object.
imshow(obj,[])  % display object
%pause
Difpat=fftshift(modulus);  % move origin
save('Difpat.mat','Difpat'); % added by Yun
save('modulus.mat','modulus'); % added by Yun
imshow(Difpat,[])  % display diffraction pattern intensity.
pause
imshow(support,[])  % display support.
%pause
phi = rand(512,512)*pi;   % first estimate of phases is random numbers.
phase = exp(i*phi);
Gkuv = modulus.*(phase);  % initial diff pattern - known amps and random phases.
gprime_kxy = ifft2(Gkuv); % first estimate of image/
gkxy = gprime_kxy.* support;  % apply known support
imshow(abs(gkxy),[]);
support2 = (support-1)*(-1);  % complement of support.

Gkuv = fft2(gkxy);
Gprime_kuv = modulus.*exp(i*angle(Gkuv));  % replace with known fourier modulus (diffracted intensities)
gprime_kxy = ifft2(Gprime_kuv);    % next extimate of image.
gk = (abs(gprime_kxy)).^2;  		% error calcn.  Note: DM program does not have square here. now same as 3D
rms1 = sum(sum(gk));				%first iteration to get rms1 (root mean square of support area after first iteration). 
gkxy = gprime_kxy.*support;			%rms1 is used in subsequent iterations as normalization factor. 
imshow(abs(gkxy),[]);
for w=1:2
 for j = 1:40                      % Do HiO.   change to desired value
   Gkuv = fft2(gkxy);
   Gprime_kuv = modulus.*exp(i*angle(Gkuv));
   gprime_kxy = ifft2(Gprime_kuv);
   gk = (abs(gprime_kxy.*support2)).^2;    % error calcn.  DM program does not have square.
   rms = sqrt(sum(sum(gk))/rms1);  % i updated this line by adding sqrt
   fprintf('HiO rms = %g.\n',rms);
   gkoutnew = gkxy.* support2 -0.8 * gprime_kxy.* support2;   % feedback factor 0.8
   gkxy = abs(gprime_kxy).*support;     % modulus constraint in real space
   gkxy = gkxy + gkoutnew;
   b1_1=abs(gkxy);          % added by Yun
   b1=['w_' int2str(w) 'j_' int2str(j)];   % added by Yun
   save(b1,'b1_1');         % added by Yun
   imshow(abs(gkxy),[]);
   pause(0.1);
 end

  for t = 1:10  % Do Error reduction 10 times.
      Gkuv = fftn(gkxy);
      Gprime_kuv = modulus.*exp(i*angle(Gkuv));  % replace with known fourier modulus (diffracted intensities)
   	gprime_kxy = ifft2(Gprime_kuv);
   	gk = (abs(gprime_kxy.*support2)).^2;  
   	rms = sqrt(sum(sum(gk))/rms1);
   	fprintf('error-reduction: rms = %g.\n',rms);
      gkxy = abs(gprime_kxy).*support;
	imshow(abs(gkxy),[]);
   end
end
imshow(abs(gkxy),[]);
