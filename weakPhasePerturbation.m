function phaseAfterPerturbation = weakPhasePerturbation(modulus,flipPercent,phase,phaseShiftAmount)
%*************************************************************************
% Convert a phse array to a new array which satisfy Friedel's law.
% You may need to use sub2ind
%*************************************************************************

%**************************************************************************
% Notation:
%     centered -> Means the direct beam/zero order Bragg spot (000) locates 
%                 at the center of 3D array ((H+1)/2,(K+1)/2,(L+1)/2)
%
%     cornered -> Means the direct beam/zero order Bragg spot (000) locates
%                 at the corner of 3D array (1,1,1)
%**************************************************************************

% Assume the modulus and angles are directly from FFT of obj, in which case
% The (000) is at the corner, rather than the center of 3D array.
modulusCentered = fftshift(modulus);
phaseCentered = fftshift(phase);

% Now search and locate the weak Bragg spots.
[modulusSort,sortIndex] = sort(modulusCentered(:));
modulusSize = size(modulusCentered);
weakIndex = sortIndex(1:floor(length(sortIndex)*flipPercent));
[weakH, weakK, weakL] = ind2sub(modulusSize,weakIndex);
weakHKL = [weakH, weakK, weakL];

% For Friedal's pair, the phi(hkl)=-phi(-h-k-l). So we only need to
% take positive numbers of either h or k or l to get unique phases.
phaseAfterPerturbationCentered = phaseCentered;
for i = 1:floor(length(weakHKL))
    % For Friedal's pair, the phi(hkl)=-phi(-h-k-l). So we only need to
    % take positive numbers of either h or k or l to get unique phases.
    if weakHKL(i,1)>=0
        weakHKL1 = weakHKL(i,:);
        weakHKL2 = modulusSize + 1 -weakHKL1;
        phaseAfterPerturbationCentered(weakHKL1) = phaseCentered(weakHKL1)+phaseShiftAmount;
        phaseAfterPerturbationCentered(weakHKL2) = phaseCentered(weakHKL2)-phaseShiftAmount;
    end
end

% Now shift the centered 3D array back to cornered 3D array which can be
% directly put in to FFT.
phaseAfterPerturbation = ifftshift(phaseAfterPerturbationCentered);