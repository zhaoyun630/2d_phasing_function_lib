function chargeSorted = sortCharge(obj)

%*************************************************************************
% This function is used to calculate the density distribution of an obj
% 
% Created by Yun Zhao on May 21th, 2015
%
% Comment: the obj can be 1D to 3D. 
% 
%*************************************************************************

% plot sorted charges over voxels
[aSize,bSize,cSize] = size(obj);
nPixels = aSize*bSize*cSize;
voxelArray = linspace(1,nPixels,nPixels);
chargeSorted = sort(obj(:));
figure
subplot(2,1,1);
plot(voxelArray,chargeSorted);
xlabel('voxels','fontsize',14);
ylabel('charges','fontsize',14);
title('sorted charge over voxel','fontsize',14);

% plot charge distribution
nInterval = 1000;
chargeCounts = zeros(nInterval+1,1);
chargeMax = chargeSorted(nPixels);
chargeMin = chargeSorted(1);
chargeInterval = (chargeMax-chargeMin)/nInterval;
for i = 1:nPixels
    countIndex = floor((chargeSorted(i)-chargeMin)/chargeInterval + 1);
    chargeCounts(countIndex) = chargeCounts(countIndex) + 1;
end

chargeArray = linspace(chargeMin,chargeMax,nInterval+1);

subplot(2,1,2);
plot(chargeArray,chargeCounts);
xlabel('charges','fontsize',14);
ylabel('counts','fontsize',14);
title('charge distribution within a unit cell','fontsize',14);
