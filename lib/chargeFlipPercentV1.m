function newDensity = chargeFlipPercentV1(objectDensity,flipPercent)

%****************************************************************
% This function is used to flip a fixed percent of density
% It sorts the charge value first, then flip lower fraction values.
% Note that charge values may be negative.
%****************************************************************

[objectDensitySorted,objectDensityIndex]=sort(objectDensity(:));
flipIndexThreshold = floor(length(objectDensityIndex)*flipPercent);
flipIndexList = objectDensityIndex(1:flipIndexThreshold);
newDensity = objectDensity;
newDensity(flipIndexList) = - objectDensity(flipIndexList);

