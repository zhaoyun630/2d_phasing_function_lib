% read the output hkl file from my python script, which extends
% the hkl to full reciprocal space by symmetry operation. The input for
% my python script is from sFall(or gen-sfs).

% first created on 3/31/2014 by Yun Zhao
% last editted on 4/3/2014 by Yun Zhao.

% This program is used to eliminate the redudant hkl rows fromt he original
% hkl file. And then fill those indices into a 3D array.

M=dlmread('3rdx_full.hkl'); % for the hkl file, you may need to delete the 
                           %  first and last line and save it as txtt.
 F = M(:,4);  % extract intensity
 F1 = 100*F/norm(F);
%  phase = M(:,5);
%  F_comp = F.*exp(1j.*phase*pi/180);
%  [L_M, W_M]=size(M);
 
 % unit cell size
 a = 57;
 b = 57;
 c = 57;

% check and delete repeated hkl line.
%  K = M(1,:);
% 
%  for i=2:20
%      ind =1;
%      for j=1:i-1
%          if M(i,1:3)== M(j,1:3)
%              ind = 0;
%          end
%      end
%      if ind ~= 0
%          K=[K;M(i,:)];
%      end
%  end
%  
%   for i=21:L_M
%      ind =1;
%      for j=i-17:i-1
%          if M(i,1:3)== M(j,1:3)
%              ind = 0;
%          end
%      end
%      if ind ~= 0
%          K=[K;M(i,:)];
%      end
%  end
 
K = M; 

dlmwrite('test.txt', K, 'precision', '%7.1f', 'newline', 'pc') % save the file
% dlmwrite('full_wo_het_3A.txt', K, 'precision', '%7.1f', 'newline', 'pc')

% put hkl list to a 3D array. 
[L_M, W_M]=size(K);
hkl = int8(K(:,1:3));
difvol_1 = zeros(a,b,c);
phase_1 = zeros(a,b,c);
angle_1 = zeros(a,b,c);
for i = 1:L_M
    h = int8(hkl(i,1)+29);
    k = int8(hkl(i,2)+29);
    l = int8(hkl(i,3)+29);
%     difvol_1(h,k,l)=K(i,4);
%     phase_1(h,k,l) = exp(1i*K(i,5)*pi/180);
    angle_1(h,k,l) = K(i,5)*pi/180;
end
% save('difvol_1.mat','difvol_1');
% save('phase_1.mat','phase_1');

save('difvol_1_full.mat','difvol_1');
save('phase_1_full.mat','phase_1');

% Now insert two zero planes along z direction
difvol_3c = zeros(a,b,3*c-2);
phase_3c = ones(a,b,3*c-2);
for i = 1:c
    num = 3*i-2; % the ith row in Difvol corresponds to 1+3(i-1) in Difvol
    difvol_3c(:,:,num)=difvol_1(:,:,i);
    phase_3c(:,:,num) = phase_1(:,:,i);
end
% save('difvol_3c.mat','difvol_3c');
% save('phase_3c.mat','phase_3c');

save('difvol_3c_full.mat','difvol_3c');
save('phase_3c_full.mat','phase_3c');

% Insert two additional zero planes in b direction

difvol_3bc = zeros(a,3*b-2,3*c-2);
phase_3bc = ones(a,3*b-2,3*c-2);

for i = 1:b
    num = 3*i-2; % the ith row in Difvol corresponds to 1+3(i-1) in Difvol
    difvol_3bc(:,num,:)=difvol_3c(:,i,:);
    phase_3bc(:,num,:) = phase_3c(:,i,:);
end

% Insert two additional zero planes in a direction

difvol_3abc = zeros(3*a-2,3*b-2,3*c-2);
phase_3abc = ones(3*a-2,3*b-2,3*c-2);

for i = 1:a
    num = 3*i-2; % the ith row in Difvol corresponds to 1+3(i-1) in Difvol
    difvol_3abc(num,:,:)=difvol_3bc(i,:,:);
    phase_3abc(num,:,:) = phase_3bc(i,:,:);
end

% save('difvol_3abc.mat','difvol_3abc');
% save('phase_3abc.mat','phase_3abc');

save('difvol_3abc_full.mat','difvol_3abc');
save('phase_3abc_full.mat','phase_3abc');



% Insert two zero planes in a,b,c direction

difvol_n = zeros(3*a-2,3*b-2,3*c-2);
phase_n = ones(3*a-2,3*b-2,3*c-2);
for i = 1:a
    for j = 1:b
        for k = 1:c
            num_i = 3*i-2; % the ith row in Difvol corresponds to 1+3(i-1) in Difvol
            num_j = 3*j-2;
            num_k = 3*k-2;
            difvol_n(num_i,num_j,num_k)=difvol_1(i,j,k);
            phase_n(num_i,num_j,num_k) = phase_1(i,j,k);
        end
    end
end
% save('difvol_n.mat','difvol_n');
% save('phase_n.mat','phase_n');

save('difvol_n_full.mat','difvol_n');
save('phase_n_full.mat','phase_n');
 







% Now creat a compace support along z direction
% support = zeros(a,b,c);
% support(:,:,25:49)=ones(a,b,c);
% save('support.mat','support');
