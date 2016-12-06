%% This part of code is used to create a phase constraint for HIO
% Created by Yun Zhao.
% Last editted by Yun Zhao on May 6th.


% load obj and calculate true phase.
obj=importdata('triple_cell.mat');
phase_true=angle(fftn(obj));

% Load the complete phase information and select planes

[a,b,c] = size(phase_true);

phase_grid_0 = zeros(a,b,c);
phase_grid_1 = ones(a,b,c);

a_c = (a+1)/2;
b_c = (b+1)/2;
c_c = (c+1)/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select a volume within a certain resolution range.

% set resolution range
resol = 0;
resol_2 = resol^2+0.001;
index_range = 2*resol + 1;

% if you want to exclude fractional Bragg orders for phase prior, use following code.
for i = 1: index_range
    for j = 1: index_range
        for k = 1: index_range
            q_h = a_c + (i-resol-1);
            q_k = b_c + (j-resol-1);
            q_l = c_c + 3*(k-resol -1);
            r_2 = (q_h-a_c)^2 + (q_k-b_c)^2 + ((q_l-c_c)/3)^2;
            if r_2 < resol_2
                phase_grid_0(q_h,q_k,q_l) = phase_grid_1(q_h,q_k,q_l);
            end
        end
    end
end


% if you want to include fractional Bragg orders for phase prior, use following code.
% a_min = a_c - resol;
% a_max = a_c + resol;
% b_min = b_c - resol;
% b_max = b_c + resol;
% c_min = c_c - 3*resol;
% c_max = c_c + 3*resol;
% 
% for i = a_min:a_max
%     for j = b_min:b_max
%         for k = c_min:c_max
%             r_2 = (i-a_c)^2 + (j-b_c)^2 + ((k-c_c)/3)^2;
%             if r_2 < resol_2
%                 phase_grid_0(i,j,k) = phase_grid_1(i,j,k);
%             end
%         end
%     end
% end



phase_true_c = fftshift(phase_true);
phase_grid_c1 = 1 - phase_grid_0;
phase_prior_c1 = phase_grid_0.*phase_true_c;
phase_prior_c = ifftshift(phase_prior_c1);
phase_grid_c = ifftshift(phase_grid_c1);

save('phase_prior.mat','phase_prior_c');
save('phase_grid.mat','phase_grid_c');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The following piece of code provide you slice selection on uncentered
% phase space from Fourier transform.

% set h = (h1,h2,...,hn) plane
% h= 20;
% for i = 1:h
%     h1 = 3*(i-1) +1;
%     phase_grid_0(:,:,h1)=phase_grid_1(:,:,h1);
% end

% phase_grid = 1 - phase_grid_0;
% phase_prior = phase_grid_0.*phase_true;
% save('phase_prior.mat','phase_prior');
% save('phase_grid.mat','phase_grid');

% for each HIO iteration, use certain know phase information.
% phase_new = phase_hio.*phase_grid + phase_prior;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% center and rearrange the phase matrix from fft.

% set h = (h1,h2,...,hn) plane
% h= 2;
% for i = 1:h
%     h1 = c_c+(i-1);
%     h2 = c_c-(i-1);
%     phase_grid_0(:,:,h1)=phase_grid_1(:,:,h1);
%     phase_grid_0(:,:,h2)=phase_grid_1(:,:,h2);
% end
% 
% 
% 
% phase_true_c = fftshift(phase_true);
% phase_grid_c1 = 1 - phase_grid_0;
% phase_prior_c1 = phase_grid_0.*phase_true_c;
% phase_prior_c = ifftshift(phase_prior_c1);
% phase_grid_c = ifftshift(phase_grid_c1);
% 
% save('phase_prior.mat','phase_prior_c');
% save('phase_grid.mat','phase_grid_c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select slices with certain angle in centered phase space.

% for plane z = 0;

% phase_grid_0(:,:,c_c) = phase_grid_1(:,:,c_c);
% 
% % for plane y = nz.(n>=1)
% ns=27;
% for h = 28-ns:4    
%     for i = 1:b
%         nz =c_c+(i-29)*3/h;
%         if mod(nz,1) == 0;
%             phase_grid_0(:,i,nz) = phase_grid_1(:,i,nz);
%         end
%     end
% 
%     % for plane x = nz.
%     for i = 1:a
%         nz = c_c+(i-29)*3/h;
%         if mod(nz,1) == 0;
%             phase_grid_0(i,:,nz) = phase_grid_1(i,:,nz);
%         end
%     end
% end
% 
% phase_true_c = fftshift(phase_true);
% phase_grid_c1 = 1 - phase_grid_0;
% phase_prior_c1 = phase_grid_0.*phase_true_c;
% phase_prior_c = ifftshift(phase_prior_c1);
% phase_grid_c = ifftshift(phase_grid_c1);
% 
% save('phase_prior.mat','phase_prior_c');
% save('phase_grid.mat','phase_grid_c');
%%




