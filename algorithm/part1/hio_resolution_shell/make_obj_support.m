% put hkl list to a 3D array.
% We cannot insert zeros along z direction ---->Yun Zhao
% Last commented by Yun Zhao at 4/17/2014

K = dlmread('full_wo_het_3A.txt');
hkl = int8(K(:,1:3));
L_M = 24526;
Difvol_1 = zeros(53,53,25);
for i = 1:L_M
    Difvol_1(hkl(i,1)+27,hkl(i,2)+27,hkl(i,3)+13)=K(i,4);
end

% Now insert two zero planes along z direction
Difvol = zeros(53,53,73);
for i = 1:25
    num = 3i-2; % the ith row in Difvol corresponds to 1+3(i-1) in Difvol
    Difvol(:,:,num)=Difvol_1(:,:,i);
end