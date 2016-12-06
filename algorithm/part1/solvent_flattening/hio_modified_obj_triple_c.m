%This piece of code is used to simulate unit cell which is applied with
%solvent flatting.


%
% UW/JCHS/YZ.  2014.  Works in Matlab R2014
% Last editted by Yun Zhao on April 17th, 2014
clear all

obj=importdata('obj_triple_c_sigma_2.mat');   
modulus1=abs(fftn(obj));    
modulus = modulus1;
support_1 = importdata('support_sigma_0_7.mat');   
[s_a,s_b,s_c] = size(support_1);

phase_true=angle(fftn(obj));
phi = rand(s_a,s_b,s_c)*2*pi;   % first estimate of phases is random numbers.
phase = exp(-1j*phi);
phase_cc_denominator=sum(sum(sum(modulus.^2)));
%phase = importdata('phase_3.mat');
Gkuv = modulus.*phase;  
gprime_kxy = ifftn(Gkuv); % first estimate of image/
gkxy = gprime_kxy.* support_1;  % apply known support
% imshow(abs(gkxy),[]);
support_2 = 1 - support_1;  % complement of support.

Gkuv = fftn(gkxy);
%Gprime_kuv = modulus.*exp(1i*angle(Gkuv));  % replace with known fourier modulus (diffracted intensities)
 Gprime_kuv = modulus.*exp(1i*(angle(Gkuv)));
gprime_kxy = ifftn(Gprime_kuv);    % next extimate of image.
gk = (abs(gprime_kxy)).^2;  		
rms1 = sum(sum(sum(gk)));			
gkxy = gprime_kxy.*support_1;			

   phase_cc_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
  phase_cc = phase_cc_numerator/phase_cc_denominator;
  rms_cc = [1 phase_cc];

for w=1:4
    for j = 1:20                      % Do HiO.   change to desired value
       Gkuv = fftn(gkxy);
       Gprime_kuv = modulus.*exp(1i*angle(Gkuv));
       gprime_kxy = ifftn(Gprime_kuv);
       gk = (abs(gprime_kxy.*support_2)).^2;    % error calcn.  DM program does not have square.
       rms = sqrt(sum(sum(sum(gk)))/rms1);  % i updated this line by adding sqrt
       phase_cc_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
%        phase_cc_numerator=sum(sum(sum(modulus.^2.*abs(cos(angle(Gkuv)-phase_true)))));
       phase_cc = phase_cc_numerator/phase_cc_denominator;
       rms_cc_info_iter = [rms phase_cc];
       rms_cc = [rms_cc; rms_cc_info_iter];
       % fprintf('HiO rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
       gkoutnew = gkxy.* support_2 -0.8 * gprime_kxy.* support_2;   % feedback factor 0.8
       gkxy = abs(gprime_kxy).*support_1;     % modulus constraint in real space
       gkxy = gkxy + gkoutnew;
       %pause(0.1);
    end

    for t = 1:20  % Do Error reduction 10 times.
        Gkuv = fftn(gkxy);
        % Gprime_kuv = modulus.*exp(1i*angle(Gkuv));  % replace with known fourier modulus (diffracted intensities)
        Gprime_kuv = modulus.*exp(1i*(angle(Gkuv)));
    	gprime_kxy = ifftn(Gprime_kuv);
    	gk = (abs(gprime_kxy.*support_2)).^2;  
        rms = sqrt(sum(sum(sum(gk)))/rms1);
        phase_cc_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
        phase_cc = phase_cc_numerator/phase_cc_denominator;
        rms_cc_info_iter = [rms phase_cc];
        rms_cc = [rms_cc; rms_cc_info_iter];
        % fprintf('error-reduction rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
    	% fprintf('error-reduction: rms = %g.\n',rms);
        gkxy = abs(gprime_kxy).*support_1;
    end
end

% obj = zeros([s_a,s_b,s_c]/3);

% Here support size is triple as the original one in each dimension.
% We just take the middle part, which is reconstructed obj.
obj_hio = gkxy(:,:,int16(s_c/3+1):int16(s_c*2/3)); 

save('obj_sf_hio.mat','obj_hio');
save('gkxy_sf.mat','gkxy');
save('rms_cc_sf.mat','rms_cc');

figure(1)
rms_cc = importdata('rms_cc_sf.mat');
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
saveas(gcf,'rms_cc_sf.png');
