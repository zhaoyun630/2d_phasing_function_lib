function newDensity = chargeFlipFixedLevel(objectDensity,flipThreshold)

%****************************************************************
% This function is used to flip a fixed percent of density
%****************************************************************


indexBelowThreshold = objectDensity<flipThreshold;
% newDensity = objectDensity;
newDensity = objectDensity.*(1-2*indexBelowThreshold);
