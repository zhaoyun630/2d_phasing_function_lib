function ccPhase = getPhaseCorrelation(phaseModel,phaseEstimated,modulusObs,modulusEst,mask)
%*************************************************************************
% Calculate Phase correlation.
%*************************************************************************


ccPhaseNumeratorMatrix = modulusEst.*modulusObs.*abs(cos(phaseEstimated-phaseModel)).*mask;
ccPhaseNumerator = sum(ccPhaseNumeratorMatrix(:));

ccPhaseDenominatorMatrix1 = (modulusObs.*mask).^2;
ccPhaseDenominatorMatrix2 = (modulusEst.*mask).^2;
ccPhaseDenominator = sqrt(sum(ccPhaseDenominatorMatrix1(:))*sum(ccPhaseDenominatorMatrix2(:)));

ccPhase = ccPhaseNumerator/ccPhaseDenominator;

% numDim = length(size(modulusObs));
% 
% switch numDim
%     case 1
%         ccPhase = sum(ccPhaseNumeratorMatrix)./sum(ccPhaseDenominatorMatrix);
%     case 2
%         ccPhase = sum(sum(ccPhaseNumeratorMatrix))./sum(sum(ccPhaseDenominatorMatrix));
%     case 3
%         ccPhase = sum(sum(sum(ccPhaseNumeratorMatrix)))./sum(sum(sum(ccPhaseDenominatorMatrix)));
%     otherwise
%         fprintf('Invalid matrix dimension. This function only works for 1~3 dimensions');
% end
