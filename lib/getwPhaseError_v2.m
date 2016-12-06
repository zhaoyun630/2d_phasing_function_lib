function weightPhaseError = getwPhaseError_v2(phaseModel,phaseEstimated,modulusObs,modulusEst,mask)
%*************************************************************************
% Calculate Weighted Phase error.
% The formular is different from v1.
% Last editted by Yun Zhao 2016.6.4
%*************************************************************************


ccPhaseNumeratorMatrix = (modulusEst+modulusObs).*acos((cos(phaseEstimated-phaseModel))).*mask;
ccPhaseNumerator = sum(ccPhaseNumeratorMatrix(:));

ccPhaseDenominatorMatrix = (modulusObs+modulusEst).*mask;
% ccPhaseDenominatorMatrix2 = (modulusEst.*mask).^2;
ccPhaseDenominator = sum(ccPhaseDenominatorMatrix(:));

weightPhaseError = ccPhaseNumerator/ccPhaseDenominator;