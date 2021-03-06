%This implementation of the Fienup HIO algorithm in 2D works for real objects. JCHS.


%
% UW/JCHS.  2003.  June 2010. Works in Matlab R2009
% Last editted by Yun Zhao, Aug 25, 2014

%2D HiO for real object (feb 03). 
% Read J.Fienup  Applied Optics 21, p.2758 (1982) and  J. Opt Soc Am. A4  p.118 (1987) to understand this program.
% Fienup hio algorithm, first written by Uwe Weierstall. TWO DIMENSIONAL.
% Save experimental modulus in DM as text file. Load it in Matlab and call it "modulus".
% e.g.: modulus = load('DMfilename'); 
% Similarly load support image and call it "support". This program has the modulus constraint in 
% real space. So it assumes a real function for the "object". Hence it would work for a weak phase object.

% obj 1x
obj =importdata('obj_2d.mat');   % This is a 512 x 512 two-dimensional array forming the test object. For electron diffraction it is the projected potential

% obj 2x
% obj = importdata('obj_2d_2x.mat');

% of a weak phase object. For xrays it is the projected charge density.
[size_a,size_b] = size(obj);
modulus=abs(fft2(obj));     %Discard phases of transform of test object. We will try to find these using the Fienup HIO algorithm.
phase_true=angle(fftn(obj));
phase_cc_denominator=sum(sum(modulus.^2));

% support 1x
%  support = importdata('support_for_acf_1x.mat');   % This is a mask, which sets to zero intensity outside the approximately known boundary of the object.
 support = importdata('support_2d.mat');
% support 2x
% support = importdata('support_for_acf_2x.mat');

% imshow(obj,[])  % display object
%pause
Difpat=fftshift(modulus);  % move origin
save('Difpat.mat','Difpat'); % added by Yun
save('modulus.mat','modulus'); % added by Yun

phi = rand(size_a,size_b)*pi;   % first estimate of phases is random numbers.
phase = exp(1i*phi);
Gkuv = modulus.*(phase);  % initial diff pattern - known amps and random phases.
gprime_kxy = ifft2(Gkuv); % first estimate of image/
gkxy = gprime_kxy.* support;  % apply known support

support2 = (support-1)*(-1);  % complement of support.

Gkuv = fft2(gkxy);
Gprime_kuv = modulus.*exp(1i*angle(Gkuv));  % replace with known fourier modulus (diffracted intensities)
gprime_kxy = ifft2(Gprime_kuv);    % next extimate of image.
gk = (abs(gprime_kxy)).^2;  		% error calcn.  Note: DM program does not have square here. now same as 3D
rms1 = sum(sum(gk));				%first iteration to get rms1 (root mean square of support area after first iteration). 
gkxy = gprime_kxy.*support;			%rms1 is used in subsequent iterations as normalization factor. 
% imshow(abs(gkxy),[]);

phase_cc_numerator=sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true)));
phase_cc = phase_cc_numerator/phase_cc_denominator;
rms_cc = [1 phase_cc];

for w=1:40
    for j = 1:30                      % Do HiO.   change to desired value
        Gkuv = fft2(gkxy);
        Gprime_kuv = modulus.*exp(1i*angle(Gkuv));
        gprime_kxy = ifft2(Gprime_kuv);
        gk = (abs(gprime_kxy.*support2)).^2;    % error calcn.  DM program does not have square.
        rms = sqrt(sum(sum(gk))/rms1);  % i updated this line by adding sqrt
%         fprintf('HiO rms = %g.\n',rms);
        gkoutnew = gkxy.* support2 -0.8 * gprime_kxy.* support2;   % feedback factor 0.8
        gkxy = abs(gprime_kxy).*support;     % modulus constraint in real space
        gkxy = gkxy + gkoutnew;
        phase_cc_numerator=sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true)));
        phase_cc = phase_cc_numerator/phase_cc_denominator;
        rms_cc = [rms_cc;rms phase_cc];
%         b1_1=abs(gkxy);          % added by Yun
%         b1=['w_' int2str(w) 'j_' int2str(j)];   % added by Yun
%         save(b1,'b1_1');         % added by Yun
%          imshow(abs(gkxy),[]);
%          pause(0.1);
    end
    
    for t = 1:10  % Do Error reduction 10 times.
        Gkuv = fftn(gkxy);
        Gprime_kuv = modulus.*exp(1i*angle(Gkuv));  % replace with known fourier modulus (diffracted intensities)
        gprime_kxy = ifft2(Gprime_kuv);
        gk = (abs(gprime_kxy.*support2)).^2;
        rms = sqrt(sum(sum(gk))/rms1);
         fprintf('error-reduction: rms = %g.\n',rms);
        gkxy = abs(gprime_kxy).*support;
        phase_cc_numerator=sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true)));
        phase_cc = phase_cc_numerator/phase_cc_denominator;
        rms_cc = [rms_cc;rms phase_cc];
%          imshow(abs(gkxy),[]);
%          pause(0.1);
    end
end
save('rms_cc_2d.mat','rms_cc');
imshow(abs(gkxy),[]);
save('gkxy.mat','gkxy');


figure(2)
rms_cc = importdata('rms_cc_2d.mat');
L = size(rms_cc);
x = linspace(1,L(1),L(1));
plot(x,rms_cc(:,1),'b-*',x,rms_cc(:,2),'r-x')
x1=xlabel('interation');
set(x1,'FontSize', 15);
set(x1,'FontWeight','bold');
y1 = ylabel('CC and rms');
set(y1,'FontSize',15);
set(y1,'FontWeight','bold')
fleg = legend('rms','CC');
set(fleg,'FontSize',15);
set(fleg,'FontWeight','bold');
saveas(gcf,'rms_cc.png');
