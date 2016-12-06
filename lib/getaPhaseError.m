function avePhaseError = getaPhaseError(phaseModel,phaseEstimated,mask)
%*************************************************************************
% Calculate Average Phase error.
%*************************************************************************


ccPhaseNumeratorMatrix = acos((cos(phaseEstimated-phaseModel))).*mask;
ccPhaseNumerator = sum(ccPhaseNumeratorMatrix(:));

ccPhaseDenominator = sum(mask(:));

avePhaseError = ccPhaseNumerator/ccPhaseDenominator;