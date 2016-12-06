function newDensity = chargeFlipPercentV2(objectDensity,flipPercent)

%****************************************************************
% This function is used to flip a fixed percent of density
% Compared with V1, it sorts by absolute value, whereas in 
% V1, it sorts by real charge value first.
%****************************************************************

objectAbsoluteDensity = abs(objectDensity);
[~,objectDensityIndex]=sort(objectAbsoluteDensity(:));
flipIndexThreshold = floor(length(objectDensityIndex)*flipPercent);
flipIndexList = objectDensityIndex(1:flipIndexThreshold);
newDensity = objectDensity;
newDensity(flipIndexList) = - objectDensity(flipIndexList);
