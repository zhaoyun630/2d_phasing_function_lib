function rFactor = getRFactor(modulusObs,modulusEstimated,mask)
%*************************************************************************
% Calculate R factor.
% Note: 
% 1, When the observed structure factor is cutted at a resolution radius 
%    and padded zeros outside of sphere to create a box. Then zero regions shouldn't 
%    be taken into account when calculating R factor between observed and estimated.
%    Becasue the observed structure factors in zeros-padded region are
%    actually unknown! If these region are taken in to account, then the R
%    factor will be very large, typically much higher than 2.
% 2, If the true structure factors are obtained by taking by FFT of model. 
%    Then it doesn't matter to impose a mask or not, in a sense that you wouldn't
%    get a ridiculous high R factors.
%*************************************************************************

% if no mask provided, set it as zero.
if mask == 0
    rFactorNumeratorMatrix=abs(modulusEstimated-modulusObs);
    rFactorDenominatorMatrix = modulusObs;

else
    rFactorNumeratorMatrix = abs(modulusEstimated-modulusObs).*mask;
    rFactorDenominatorMatrix = modulusObs.*mask;
end


numDim = length(size(modulusObs));
% fprintf('numDim is %d\n', numDim);

switch numDim
    case 1
        rFactor = sum(rFactorNumeratorMatrix)./sum(rFactorDenominatorMatrix);
    case 2
        rFactor = sum(sum(rFactorNumeratorMatrix))./sum(sum(rFactorDenominatorMatrix));
    case 3
        rFactor = sum(sum(sum(rFactorNumeratorMatrix)))./sum(sum(sum(rFactorDenominatorMatrix)));
    otherwise
        fprintf('Invalid matrix dimension. This function only works for 1~3 dimensions');
end
