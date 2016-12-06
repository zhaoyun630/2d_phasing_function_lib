%This implementation of the Fienup HIO algorithm in 3D works for real objects. JCHS.


%
% UW/JCHS/YZ.  2014.  Works in Matlab R2014
% Last editted by Yun Zhao on April 21th, 2014
% Compared with version 2, structure at each iteration is recorded in
% current version.

clear all

obj=importdata('triple_cell.mat');
phase_true=angle(fftn(obj));

% Load the complete phase information and select planes

[a,b,c] = size(obj);
a_c = (a+1)/2;
b_c = (b+1)/2;
c_c = (c+1)/2;
n=2;
a_r = (a_c-n):(a_c+n);
b_r = (b_c-n):(b_c+n);
c_r = (c_c-n):(c_c+n);



% select slices with certain angle in centered phase space.

% for plane z = 0;

% phase_grid_0(:,:,c_c) = phase_grid_1(:,:,c_c);
% 
% % for plane y = nz.(n>=1)

for h = 1:28
    phase_grid_0 = zeros(a,b,c);
    phase_grid_1 = ones(a,b,c);
    for i = 1:b
        nz =c_c+(i-29)*3/h;
        if mod(nz,1) == 0;
            phase_grid_0(:,i,nz) = phase_grid_1(:,i,nz);
        end
    end
    
    % % for plane x = nz.
    %     for i = 1:a
    %         nz = c_c+(i-29)*3/h;
    %         if mod(nz,1) == 0;
    %             phase_grid_0(i,:,nz) = phase_grid_1(i,:,nz);
    %         end
    %     end
    
    
    phase_true_c = fftshift(phase_true);
    phase_grid_c1 = 1 - phase_grid_0;
    phase_prior_c1 = phase_grid_0.*phase_true_c;
    phase_prior_c = ifftshift(phase_prior_c1);
    phase_grid_c = ifftshift(phase_grid_c1);
    
    save('phase_prior.mat','phase_prior_c');
    save('phase_grid.mat','phase_grid_c');
    
    
    % obj=importdata('triple_cell.mat');
    modulus1=abs(fftn(obj));     %Discard phases of transform of test object. We will try to find these using the Fienup HIO algorithm.
    %modulus = poissrnd(modulus1);
    
    
    modulus = modulus1;
    support = importdata('support_3d.mat');   % This is a mask, which sets to zero intensity outside the approximately known boundary of the object.
    phase_prior = importdata('phase_prior.mat');
    phase_grid = importdata('phase_grid.mat');
    
    %phase_true=angle(fftn(obj));
    phi_1 = rand(a,b,c)*2*pi;   % first estimate of phases is random numbers.
    phi = phi_1.*phase_grid + phase_prior;
    %   phi = phi_1;
    % phi = phase_true;
    %phi = poissrnd(phase_true);
    %phi = zeros(53,53,73);
    phase = exp(1i*phi);
    % phase_cc_denominator=sum(sum(sum(modulus.^2)));
    phase_cc_denominator=sum(sum(sum(modulus.^2)))-modulus(a_c,b_c,c_c)^2;
    %  phase_cc_denominator=sum(sum(sum(modulus.^2)))-sum(sum(sum(modulus(a_r,b_r,c_r).^2)));
    %phase = importdata('phase_3.mat');
    Gkuv = modulus.*phase;  % initial diff pattern - known amps and random phases.
    phase_w = modulus.^2.*cos(angle(Gkuv)-phase_true);
    phase_cc_numerator=sum(sum(sum(phase_w)))-phase_w(a_c,b_c,c_c);
    %   phase_cc_numerator=sum(sum(sum(phase_w)))-sum(sum(sum(phase_w(a_c,b_c,c_c))));
    %   phase_cc_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
    phase_cc = phase_cc_numerator/phase_cc_denominator;
    rms_cc = [1 phase_cc];
    
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
    %     phase_cc_numerator=sum(sum(sum(modulus.^2.*abs(cos(angle(Gkuv)-phase_true)))));
    phase_w = modulus.^2.*cos(angle(Gkuv)-phase_true);
    phase_cc_numerator=sum(sum(sum(phase_w)))-phase_w(a_c,b_c,c_c);
    % phase_cc_numerator=sum(sum(sum(phase_w)))-sum(sum(sum(phase_w(a_c,b_c,c_c))));
    % phase_cc_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
    phase_cc = phase_cc_numerator/phase_cc_denominator;
    rms_cc = [rms_cc; 1 phase_cc];
    % fprintf('HiO rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
    for w=1:4
        for j = 1:20                      % Do HiO.   change to desired value
            Gkuv = fftn(gkxy);
            % Gprime_kuv = modulus.*exp(1i*angle(Gkuv));
            Gprime_kuv = modulus.*exp(1i*(angle(Gkuv).*phase_grid +phase_prior));
            gprime_kxy = ifftn(Gprime_kuv);
            gk = (abs(gprime_kxy.*support2)).^2;    % error calcn.  DM program does not have square.
            rms = sqrt(sum(sum(sum(gk)))/rms1);  % i updated this line by adding sqrt
            % phase_cc_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
            %   phase_cc_numerator=sum(sum(sum(modulus.^2.*abs(cos(angle(Gkuv)-phase_true)))));
            phase_w = modulus.^2.*cos(angle(Gkuv)-phase_true);
            phase_cc_numerator=sum(sum(sum(phase_w)))-phase_w(a_c,b_c,c_c);
            % phase_cc_numerator=sum(sum(sum(phase_w)))-sum(sum(sum(phase_w(a_c,b_c,c_c))));
            % phase_cc_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
            phase_cc = phase_cc_numerator/phase_cc_denominator;
            rms_cc_info_iter = [rms phase_cc];
            rms_cc = [rms_cc; rms_cc_info_iter];
            % fprintf('HiO rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
            gkoutnew = gkxy.* support2 -0.8 * gprime_kxy.* support2;   % feedback factor 0.8
            gkxy = abs(gprime_kxy).*support;     % modulus constraint in real space
            
            %       b1_1=abs(gkxy);          % added by Yun
            %    b1=['w_' int2str(w) 'j_' int2str(j)];   % added by Yun
            %    save(b1,'b1_1');         % added by Yun
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
            %phase_cc_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
            % phase_cc_numerator=sum(sum(sum(modulus.^2.*abs(cos(angle(Gkuv)-phase_true)))));
            %     phase_w = modulus.^2.*cos(angle(Gkuv)-phase_true);
            %   phase_cc_numerator=sum(sum(sum(phase_w)))-phase_w(a_c,b_c,c_c);
            % phase_cc_numerator=sum(sum(sum(phase_w)))-sum(sum(sum(phase_w(a_c,b_c,c_c))));
            phase_cc_numerator=sum(sum(sum(modulus.^2.*cos(angle(Gkuv)-phase_true))));
            phase_cc = phase_cc_numerator/phase_cc_denominator;
            
            rms_cc_info_iter = [rms phase_cc];
            rms_cc = [rms_cc; rms_cc_info_iter];
            % fprintf('error-reduction rms = %g; \t phase_error_perc = %g.\n',rms, phase_error_perc);
            %fprintf('error-reduction: rms = %g.\n',rms);
            gkxy = abs(gprime_kxy).*support;
            %          b1_1=abs(gkxy);          % added by Yun
            %    b1=['w_' int2str(w) 't_' int2str(t)];   % added by Yun
            %    save(b1,'b1_1');         % added by Yun
            
        end
    end
    
    rms_cc_h = ['rms_cc_' int2str(h) '.mat'];
    gkxy_h = ['gkxy_' int2str(h) '.mat'];
    save(gkxy_h,'gkxy');
    save(rms_cc_h,'rms_cc');
end

% plot rms and CC
figure(1)
rms_cc_h = ['rms_cc_' int2str(h) '.mat'];
rms_cc = importdata(rms_cc_h);
L = size(rms_cc);
x = linspace(1,L(1),L(1));
plot(x,rms_cc(:,1),'b',x,rms_cc(:,2),'r')
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

