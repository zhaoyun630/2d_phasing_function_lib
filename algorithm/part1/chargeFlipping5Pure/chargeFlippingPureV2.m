%######################################################################
%       Ab initio phasing with pure charge flipping algorithm
%######################################################################
%
% First created by Yun Zhao on Nov 12th, 2014 at Xiaotao, China
% Last editted by Yun Zhao on May 15th, 2015
% 
%*******************************************************************
% About Version 2:
%       Version 1 is successful. It has preliminary success.
%       In version2, we refined the algorithm to obtain better
%       reconstruction. Especially, we modify the flipFrac during
%       iteration. The main index is Rfactor
%
%*************************************************************************
% 1) This piece of code exclusively use charge flipping
% 2) Notation:
%     centered -> Means the direct beam/zero order Bragg spot (000) locates 
%                 at the center of 3D array ((H+1)/2,(K+1)/2,(L+1)/2)
%
%     cornered -> Means the direct beam/zero order Bragg spot (000) locates
%                 at the corner of 3D array (1,1,1)
%**************************************************************************

% Whenever you start a new run, set a new runID first.
runID = 8;

% Using Phase Perturbation method? Set 1 if you want to use it.
phasePerturbation = 0;

% Flip at a fixed value(set 0) or fixed percentage(set 1).
flipMethod = 1;

% set some global variables.
modelFile = 'unit_cell_shift_1A.mat';
reflectionlist = 'recLattice1A.mat';
% for spinel: 'HKLspinelN25.mat';
% for streptavidin: 'recLattice1A.mat';
flipChargeLevel = 0.8;
flipFrac = 0.8;
weakAmpPercent = 0.1;
weakPhaseShift = pi/2;
iterNumTotal = 5000;
runNotes = ['The object is a streptavidin crystal. '...
            'Flip at fixed percentage level. '...
            'Weak phase perturbation is defaulted. '...
            'Hope it is successful'];

runConditionFile = strcat('runID_',int2str(runID),'_runCondition.mat');
save(runConditionFile,'runID','phasePerturbation','flipMethod','reflectionlist','flipChargeLevel','flipFrac','weakAmpPercent','weakPhaseShift','iterNumTotal','runNotes');

% Load 3D reflection array.
% Note that reflection list is actually complex structure factors.
objModel = importdata(modelFile); % it is only used to search origin.
FComplexCenteredObserved = importdata(reflectionlist);
FModulusCenteredObserved = abs(FComplexCenteredObserved);
FModulusCorneredObserved = ifftshift(FModulusCenteredObserved);

% Create a mask for reflection list. We only take Bragg spots which are
% observed in account when calculating R factors and other metrics.
FMaskCentered = FModulusCenteredObserved > 0.01;
FMaskCornered = ifftshift(FMaskCentered);

% Real object and its size
obj = real(ifft(ifftshift(FComplexCenteredObserved)));
[sa,sb,sc] = size(obj);

% Create an array to store R factors.
RFactor = zeros(iterNumTotal,1);

% Number of voxels to be flipped.
flipIndex = floor(flipFrac*sa*sb*sc);

% Creat random phase with centrosymmetry.
phase0 = friedelPhaseRandom([sa,sb,sc],'cornered');

% First estimation of our complex structure factors from random phases
FComp0 = FModulusCorneredObserved.*exp(1i*phase0);

% First estimation of our object from random phases.
objK = real(ifftn(FComp0));

flipFracIter = flipFrac;
% Start Charge flipping iteration.
for iter = 1: iterNumTotal
    % Charge flipping in real space
    if flipMethod == 1
        objKFlip = chargeFlipPercent(objK,flipFracIter);
    elseif flipMethod == 0;
        objKFlip = chargeFlipFixedLevel(objK,flipChargeLevel);
    else
        fprintf('Please choose a valid flip method. Flipping fixed charge(set 0) or fixed percentage (set 1)\n');
    end
    
    % Constraints in Fourier Space
    FComp1 =fftn(objKFlip);
    FModulusIter= abs(FComp1);
    phaseIter = angle(FComp1);
    
    RFactor(iter) = getRFactor(FModulusCorneredObserved,FModulusIter,FMaskCornered);   
    FModulus = modulusConstraint(FModulusCorneredObserved,FModulusIter,FMaskCornered);
    
    if RFactor(iter) >= 0.1 && RFactor(iter)<= 0.35     
        flipFracIter = 0.999*flipFracIter;
    end
%     RFactor(iter) = getRFactor(FModulusCorneredObserved,FModulusIter,0);   
%     FModulus = modulusConstraint(FModulusCorneredObserved,FModulusIter,0); 
    if mod(iter,100) == 0
        fprintf('%dth iteration. R factor %f \n',iter,RFactor(iter));
    end

    % Do weak phase perturbation?
    if phasePerturbation == 1
        phasePerturbed = weakPhasePerturbation(FModulusCorneredObserved,flipFracIter,phaseIter,weakPhaseShift);
        phaseIter = phasePerturbed;
%     else
%         fprintf('set a valid phasePerturbation method \n');
    end
    
    % Combine measured amplitudes with estimated phase
    F_comp_2 = FModulus.*exp(1i*phaseIter);
    objK = real(ifftn(F_comp_2));
    
end

% Construct a object to save obj estimate and relavant parameters.
originShift = getOrigin(objK,objModel);
objCentered = circshift(objK,originShift);

objName = strcat('runID_',int2str(runID),'_obj.mat');
metricName = strcat('runID_',int2str(runID),'_metric.mat');
save(objName,'objCentered');
save(metricName,'RFactor');

% Plot Figures
RPhaseCCPlot(RFactor,runID);

fprintf('Now show the 3D view of reconstructed object: \n')
% Show 3D view.
show3DObjv2(objK,0.998,1,0.9)



