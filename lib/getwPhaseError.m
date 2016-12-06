function weightPhaseError = getwPhaseError(phaseModel,phaseEstimated,modulusObs,modulusEst,mask)
%*************************************************************************
% Calculate Weighted Phase error.
%*************************************************************************


ccPhaseNumeratorMatrix = modulusEst.*modulusObs.*acos((cos(phaseEstimated-phaseModel))).*mask;
ccPhaseNumerator = sum(ccPhaseNumeratorMatrix(:));

ccPhaseDenominatorMatrix1 = (modulusObs.*mask).^2;
ccPhaseDenominatorMatrix2 = (modulusEst.*mask).^2;
ccPhaseDenominator = sqrt(sum(ccPhaseDenominatorMatrix1(:))*sum(ccPhaseDenominatorMatrix2(:)));

weightPhaseError = ccPhaseNumerator/ccPhaseDenominator;