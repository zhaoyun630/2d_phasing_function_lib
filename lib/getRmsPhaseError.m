function rmsPhaseError = getRmsPhaseError(phaseModel,phaseEstimated,modulusObs,modulusEst,mask)
%*************************************************************************
% Calculate rms Phase error (Frank,1996). John used this metric in (Spence et al,2003).
%*************************************************************************


ccPhaseNumeratorMatrix = (modulusEst+modulusObs).*(phaseEstimated-phaseModel).^2.*mask;
ccPhaseNumerator = sum(ccPhaseNumeratorMatrix(:));

ccPhaseDenominatorMatrix1 = modulusObs.*mask;
ccPhaseDenominatorMatrix2 = modulusEst.*mask;
ccPhaseDenominator = sum(ccPhaseDenominatorMatrix1(:)) + sum(ccPhaseDenominatorMatrix2(:));

rmsPhaseError = sqrt(ccPhaseNumerator/ccPhaseDenominator);