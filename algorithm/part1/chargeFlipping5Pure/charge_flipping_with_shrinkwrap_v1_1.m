% The algorithm implemented in this code:
% In each cycle:
%   1) Charge flipping
%   2) several HIO+ER cycle
%       2.1) HIO
%       2.2) ER
%       2.3) Shrinkwrap

% First created by Yun Zhao on Mar 4, 2015
% Last editted by Yun Zhao on Mar 4, 2015


% offcentered means F(000) is located at the corner of 3D array, not center
F_modulus_offcentered_observed = importdata('rec_tri_lat_mag_1A.mat');
support_1 = importdata('support_triple_1A.mat');
support_2 = 1-support_1;  % complement of support.

% Set parameters for Guassian filtering during shrinkwrap step
guassian_sigma = 2.5; % The initial sigma is 3 pixels
guassian_mask_threshold = 0.05; % The cutoff is set at 20% of its maximum value

% set some global variables.
solvent_level = 0.7;
flip_frac = 79;
weak_amp_percent = 0.001;
R_factor = [];
F_000 = [];
[s_a,s_b,s_c] = size(F_modulus_offcentered_observed);
flip_index = floor(flip_frac/100*s_a*s_b*(s_c+2)/3);

% Iteration parameters
data_resolution = '1A';
run_number = 'v1';
iter_big_cycle = 3;
iter_charge_flipping = 200;
iter_HIO_ER_cycle = 3;
iter_HIO = 40;
iter_ER = 20;

% name appear in string. It will be used to label output file name.
num2str1 = int2str(flip_frac);
num2str2 = int2str(iter_big_cycle);
num2str3 = int2str(iter_charge_flipping);
num2str4 = int2str(iter_HIO_ER_cycle);
num2str5 = int2str(iter_HIO);
num2str6 = int2str(iter_ER);

parameter_name = strcat(data_resolution,'_',num2str1,'_',num2str2,'_',num2str3,'_',num2str4,'_',num2str5,'_',num2str6,'_',run_number);



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
modulus_centered_0 = fftshift(F_modulus_offcentered_observed);
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
obj_estimate = real(ifftn(ifftshift(F_comp_0)));

% % Start Charge flipping iteration
% obj_1d_array = reshape(obj_k,[1,s_a*s_b*s_c]);
% rho_sigma = std(obj_1d_array);
% rho_mean = mean(obj_1d_array);
% rho_max = max(obj_1d_array);
% 
% 
% % rho_cut_off = solvent_level*rho_sigma;
% rho_sort = sort(obj_1d_array);
% % rho_cut_off = rho_sort(flip_index);
% % rho_cut_off = 0.01;

for iter_cycle = 1:iter_big_cycle
    % Start with charge flipping to find the boundary in a-b plane
    for iter_c = 1: iter_charge_flipping
        % Charge flipping in real space
        obj_k_middle = obj_estimate(:,:,(s_c+2)/3:(2*s_c+1)/3); % Actually I apply and compact support along z here.
        obj_1d_array = reshape(obj_k_middle,[1,s_a*s_b*(s_c+2)/3]);
        rho_sigma = std(obj_1d_array);
        rho_mean = mean(obj_1d_array);
        rho_max = max(obj_1d_array);
        rho_sort = sort(obj_1d_array);
        rho_cut_off = rho_sort(flip_index);
        %     rho_cut_off = solvent_level*rho_sigma;
        obj_k_middle_flip= obj_k_middle;
        for i = 1:s_a
            for j = 1:s_b
                for k = 1:(s_c+2)/3
                    if abs(obj_k_middle(i,j,k))<rho_cut_off
                        obj_k_middle_flip(i,j,k) = -obj_k_middle(i,j,k);
                    end
                end
            end
        end
        
        % Now pad zeros above and below to get a triple cell.
        obj_k_1 = zeros(s_a,s_b,s_c);
        obj_k_1(:,:,(s_c+2)/3:(2*s_c+1)/3) = obj_k_middle_flip;
        
        % Constraints in Fourier Space
        F_comp_1 =fftn(obj_k_1);
        F_modulus= abs(F_comp_1);
        F_modulus_constraint = F_modulus_offcentered_observed;
        F_modulus_constraint(1,1,1)=F_modulus(1,1,1);
        r_i = sum(sum(sum(abs(abs(F_modulus)-abs(F_modulus_constraint)))))/(sum(sum(sum(abs(F_modulus_constraint))))-abs(F_modulus_constraint(1,1,1)));
        R_factor = [R_factor;r_i];
        F_000 = [F_000; F_modulus(1,1,1)];
        phase = angle(F_comp_1);
        
        
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
        obj_estimate = real(ifftn(F_comp_2));
        
    end
    
    obj_estimate = obj_estimate.*support_1;
    
    
    
    % Now save the intermediate structure after charge flipping
    obj_estimate_one_unit_inter = obj_estimate(:,:,(s_c+2)/3:(2*s_c+1)/3);
    obj_intermediate_name = strcat('obj_estimate_',int2str(iter_cycle),'_intermediate_',parameter_name,'.mat');
    save(obj_intermediate_name,'obj_estimate_one_unit_inter');
    
 
    
    % Now apply HIO iteration.
    for w=1:iter_HIO_ER_cycle
        for j = 1:iter_HIO                      % Do HiO.   change to desired value
            F_comp_1 = fftn(obj_estimate);
            F_comp_2 = F_modulus_offcentered_observed.*exp(1i*angle(F_comp_1));
            obj_estimate_prime = real(ifftn(F_comp_2));
            obj_outside_support = (real(obj_estimate_prime.*support_2)).^2;    % error calcn.  DM program does not have square.
            rms = sqrt(sum(sum(sum(obj_outside_support)))/sum(sum(sum(abs(obj_estimate_prime)))));  % i updated this line by adding sqrt
            F_modulus = abs(F_comp_1);
            F_000 = [F_000; F_modulus(1,1,1)];
            F_modulus_constraint = F_modulus_offcentered_observed;
            F_modulus_constraint(1,1,1)=F_modulus(1,1,1);
%             F_modulus(1,1,1) = F_modulus_constraint(1,1,1);
            r_i = sum(sum(sum(abs(abs(F_modulus)-abs(F_modulus_constraint)))))/(sum(sum(sum(abs(F_modulus_constraint))))-abs(F_modulus_constraint(1,1,1)));
            R_factor = [R_factor;r_i];
            obj_out_new = obj_estimate.* support_2 -0.8 * obj_estimate_prime.* support_2;   % feedback factor 0.8
            obj_estimate_1 = real(obj_estimate_prime).*support_1;     % modulus constraint in real space
            obj_estimate = real(obj_estimate_1 + obj_out_new);
        end
        
        for t = 1:iter_ER  % Do Error reduction 10 times.
            F_comp_1 = fftn(obj_estimate);
            F_comp_2 = F_modulus_offcentered_observed.*exp(1i*(angle(F_comp_1)));
            obj_estimate_prime = real(ifftn(F_comp_2));
            obj_outside_support = (abs(obj_estimate_prime.*support_2)).^2;
            rms = sqrt(sum(sum(sum(obj_outside_support)))/sum(sum(sum(abs(obj_estimate_prime)))));
            F_modulus = abs(F_comp_1);
            F_000 = [F_000; F_modulus(1,1,1)];
            F_modulus_constraint = F_modulus_offcentered_observed;
            F_modulus_constraint(1,1,1)=F_modulus(1,1,1);
%             F_modulus(1,1,1) = F_modulus_constraint(1,1,1);
            r_i = sum(sum(sum(abs(abs(F_modulus)-abs(F_modulus_constraint)))))/(sum(sum(sum(abs(F_modulus_constraint))))-abs(F_modulus_constraint(1,1,1)));
            R_factor = [R_factor;r_i];
            obj_estimate = real(obj_estimate_prime).*support_1;
        end
        
        
            % Now perform shrinkwrap.
            % Firstly, convolute the density with a guassian.
            % Set a guassian filter
            if guassian_sigma < 1.5
                guassian_sigma = guassian_sigma*(101-((iter_cycle-1)*iter_HIO_ER_cycle+w))/100;
            else
                guassian_sigma = 1.5;
            end
            guassian_filter_size = guassian_sigma;
            mask_coord= -guassian_filter_size:1:guassian_filter_size;
            [gmask_x,gmask_y,gmask_z] = meshgrid(mask_coord,mask_coord,mask_coord);
            Guassian_filter_prime = normpdf(gmask_x,0,guassian_sigma).*normpdf(gmask_y,0,guassian_sigma).*normpdf(gmask_z,0,guassian_sigma);
            Guassian_filter = Guassian_filter_prime/sum(sum(sum(Guassian_filter_prime)));
            
            % Now convolute the object density with a guassian filter
            obj_estimate_conv_guass = imfilter(obj_estimate,Guassian_filter);
            
            % Calculate some statistic numbers for support update
            obj_estimate_1d_convol_guass = reshape(obj_estimate_conv_guass,[1,s_a*s_b*s_c]);
            %     obj_estimate_sigma = std(obj_estimate_1d_convol_guass);
            %     obj_estimate_mean = mean(obj_estimate_1d_convol_guass);
            obj_estimate_cut_off = guassian_mask_threshold*max(obj_estimate_1d_convol_guass);
            
            % Apply the cut off value to shrink support.
            support_1d = zeros(1,s_a*s_b*s_c);
            for i = 1:s_a*s_b*s_c
                if obj_estimate_1d_convol_guass(i)>=obj_estimate_cut_off
                    support_1d(i) = 1;
                end
            end
            support_1 = reshape(support_1d,[s_a,s_b,s_c]);
            %     support_1 = reshape(support_1d,[s_a,s_b,s_c]).* support_0;
            support_2 = 1- support_1; % complimentary support
    end   
end


% Note Here I name the file in following way
% The first   : solvent fraction
% The second  : # of big cf+HIO cycle
% The third   : # of charge flipping cycle
% The forth   : # of big HIO cycle
% The fifth   : # of HIO cycle
% The sixth   : # of ER cycle
% The seventh : # of run

% The structure reconstructed locates at the middle of triple cell.
obj_estimate_one_unit = obj_estimate(:,:,(s_c+2)/3:(2*s_c+1)/3);

obj_file_name = strcat('obj_estimate_',parameter_name,'.mat');
r_factor_file_name = strcat('r_factor_',parameter_name,'.mat');
F_000_file_name = strcat('F_000_',parameter_name,'.mat');
% save('obj_estimate_1A_79_3_100_3_40_10_1.mat','obj_estimate_one_unit');
% save('r_factor_1A_79_3_100_3_40_10_1.mat','R_factor');
% save('F_000_1A_79_3_100_3_40_10_1.mat','F_000');
save(obj_file_name,'obj_estimate_one_unit');
save(r_factor_file_name,'R_factor');
save(F_000_file_name,'F_000');


% R_factor = importdata('r_factor_file_name');
% F_000 = importdata('F_000_file_name');
iter_num_total = length(R_factor);
iter_num = linspace(1,iter_num_total,iter_num_total);
plt1 = figure(1);
plot(iter_num,R_factor,'-b.');
% lg1_text1 = strcat('total iteration: ',int2str(iter_num_total)); 
% lg1_text2 = strcat('weak reflection percentage:' , int2str(weak_amp_percent));
% legend(lg1_text1,lg1_text2);
% text(700,1.54,lg1_text1,'FontSize', 15,'FontWeight','bold');
% text(700,1.53,lg1_text2,'FontSize', 15,'FontWeight','bold');
xlabel('iteration','FontSize', 15,'FontWeight','bold');
ylabel('R factor','FontSize', 15,'FontWeight','bold');
% rectangle('Position',[690,1.525,250,0.022]);
r_factor_fig_name = strcat('r_factor_',parameter_name,'.tif');
saveas(plt1,r_factor_fig_name);

plt2 = figure(2);
plot(iter_num,F_000,'-r.');
% lg1_text1 = strcat('total iteration: ',int2str(iter_num_total)); 
% lg1_text2 = strcat('weak reflection percentage:' , int2str(weak_amp_percent));
% legend(lg1_text1,lg1_text2);
% text(700,1.54,lg1_text1,'FontSize', 15,'FontWeight','bold');
% text(700,1.53,lg1_text2,'FontSize', 15,'FontWeight','bold');
xlabel('iteration','FontSize', 15,'FontWeight','bold');
ylabel('F(000)','FontSize', 15,'FontWeight','bold');
% rectangle('Position',[690,1.525,250,0.022]);
F_000_fig_name = strcat('F_000',parameter_name,'.tif');
saveas(plt2,F_000_fig_name);




