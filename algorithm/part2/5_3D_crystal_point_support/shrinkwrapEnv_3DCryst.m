% ######################################################################
%  HIO alogrithm given with a molecular envelope for 3D crystal
% ######################################################################
%% Version
% First Created by Yun Zhao on Feb 13rd, 2015
% Last editted by Yun Zhao on Jun 2nd, 2016
% 
% Compared with shrinkwrapEnv_v1_*.m or shrinkwrapEnv.m, This version is
% used for 3D crystal dataset. 
% *************************************************************************
%% Notation:
%     centered -> Means the direct beam/zero order Bragg spot (000) locates 
%                 at the center of 3D array ((H+1)/2,(K+1)/2,(L+1)/2)
%
%     cornered -> Means the direct beam/zero order Bragg spot (000) locates
%                 at the corner of 3D array (1,1,1)
%*********************************************************

%% Input information
% Whenever you start a new run, set a new runID first.
runID = 9;

% Input ID
inputID = 9;

% default Shrinkwrap. Set 1 if you want to use it.
shrinkwrap = 0;

% The hole constraint is always imposed? Set 1 if you want it.
constantSupportConstraint = 0;

% Impose symmetry constraint: set 1;
symConstraint = 0;

% Get reflection list from model with real numbers(set 1) 
% or abs values(2) or experimental value (set 0)
reflectionListSource = 1;

% set some global variables.
% modelFile = 'inputID2_newObj.mat';
modelFile = strcat('inputID',int2str(inputID),'_newObj.mat');

% reflectionlist = 'HKLTriplecellCentered3A_v2.mat';
% reflectionlist = 'inputID2_fcomplex.mat';
reflectionlist = strcat('inputID',int2str(inputID),'_fcomplex.mat');

% supportFile = 'hole50Perc3A.mat';
% supportFile ='inputID2_triSupport.mat';
supportFile = strcat('inputID',int2str(inputID),'_Support.mat');

betaHIO = 0.8; % feedback parameter in HIO

% solventPercent = 0.7;
% solventDensity = 0;
guassianSigma = 2.5; % The initial sigma is 3 pixels
guassianSigmaMin = 1.5;
guassianMaskThreshold = 0.05; % The cutoff is set at 5% of its maximum value

iterComb = 50;
iterHIO = 30;
iterER = 10;
runNotes = ['The object is a triple cell. '...
            'The protein pdb is 3rdx'];

runConditionFile = strcat('runID_',int2str(runID),'_runCondition.mat');

save(runConditionFile,'runID','reflectionlist','supportFile','betaHIO','guassianSigma','guassianSigmaMin','guassianMaskThreshold','iterComb','iterHIO','iterER','runNotes');

% Load model
model = importdata(modelFile);
model_s_abs = abs(fftn(model));
model_s_angle = angle(fftn(model));

% Load 3D reflection array.
% Note that reflection list is actually complex structure factors.
FComplexCorneredObserved = importdata(reflectionlist);

% get modules and phases from model
FModulusCorneredModel = abs(FComplexCorneredObserved);
FPhaseCorneredModel = angle(FComplexCorneredObserved);
% FModulusCorneredModel = ifftshift(FModulusCenteredModel);
% FPhaseCorneredModel=ifftshift(FPhaseCenteredModel);

% HKLmask = FModulusCorneredModel > 0.01;
HKLmask = model_s_abs > 0.01;

% This is initial support, which sets to zero intensity outside the approximately known boundary of the object.
support0 = importdata(supportFile);   
support1 = support0;
[sizeA,sizeB,sizeC] = size(support1);

% Create an array to record rms, cc, Rfactor over iteration.
figureOfMerit = zeros(iterComb*(iterHIO+iterER),5);


%% Initialization
% The code starts from here.
% first estimate of phases is random numbers.
FPhaseInitial = rand(sizeA,sizeB,sizeC)*2*pi;   
FPhaseExpTerm = exp(1j*FPhaseInitial);
FComplexCorneredEstimatePrime = FModulusCorneredModel.*FPhaseExpTerm; % Complex structure factor amplitude

% The following line is used to calculate CC value. The donominator value is same for all iterations.
phaseCCDenominator=sum(sum(sum(FModulusCorneredModel.^2)));
support2 = 1-support1;  % complement of support.

% first estimate of object
objEstimatePrime = abs(ifftn(FComplexCorneredEstimatePrime)); 

% calculate first estimate of figure of merits

% Initial estimate of model in unit cell, |F| and phase
% objIter = objEstimatePrime(:,:,int32(sizeC/3+1):int32(sizeC*2/3));
objIter = objEstimatePrime;
objIter_s_abs = abs(fftn(objIter));
objIter_s_phase = angle(fftn(objIter));

% for unit cell
% Calculate rms
iterIndex = 1;
obj_outside_support = objEstimatePrime.*support2;    % error calcn.  DM program does not have square.
rms = sqrt(sum(sum(sum(obj_outside_support.^2)))/sum(sum(sum(abs(objEstimatePrime.^2)))));  % i updated this line by adding sqrt
figureOfMerit(iterIndex,1) = rms;
% Calculate cc
figureOfMerit(iterIndex,2) = getPhaseCorrelation(model_s_angle,objIter_s_phase,model_s_abs,objIter_s_abs,HKLmask);
% Calculate R_factor
figureOfMerit(iterIndex,3) = getRFactor(model_s_abs,objIter_s_abs,HKLmask);
% Calculate average phase_error
figureOfMerit(iterIndex,4) = getaPhaseError(model_s_angle,objIter_s_phase,HKLmask);
% Calculate weighted phase_error
figureOfMerit(iterIndex,5) = getwPhaseError(model_s_angle,objIter_s_phase,model_s_abs,objIter_s_abs,HKLmask);


% First estimate of object after applying support
objEstimate = abs(objEstimatePrime.* support1);  % apply known support


%% Start HIO and Shrinkwrap
% objEstimate = model;

for i=1:iterComb
    % Do HiO.   change to desired value
    for j = 1:iterHIO
        iterIndex = (i-1)*(iterHIO+iterER) + j+1;
        FComplexCorneredEstimatePrime = fftn(objEstimate);
%         FModulusCorneredEstimate = abs(FComplexCorneredEstimatePrime);
        FPhaseCorneredEstimate = angle(FComplexCorneredEstimatePrime);
        
        % Current estimate of model in unit cell, |F| and phase
%         objIter = objEstimate(:,:,int32(sizeC/3+1):int32(sizeC*2/3));
        objIter = objEstimate;
        objIter_s_abs = abs(fftn(objIter));
        objIter_s_phase = angle(fftn(objIter));
        
        % The modulus constraint is implemented in the following line
        FComplexCorneredEstimated = FModulusCorneredModel.*exp(1i*FPhaseCorneredEstimate);
        
        % Inverse Fourier Transform gives object.
        objEstimatePrime = abs(ifftn(FComplexCorneredEstimated));
        
        % Impose support from real space
        objOutsideSupport = objEstimate.* support2 - betaHIO * objEstimatePrime.* support2;   % feedback factor 0.8
        objEstimate1 = (objEstimatePrime).*support1;     % modulus constraint in abs space
        objEstimate = abs(objEstimate1) + abs(objOutsideSupport); % new object estimation for next iteration.
        
        
        % Calculate rms - figure of merits
        obj_outside_support = objEstimatePrime.*support2;    % error calcn.  DM program does not have square.
        rms = sqrt(sum(sum(sum(obj_outside_support.^2)))/sum(sum(sum(abs(objEstimatePrime.^2)))));  % i updated this line by adding sqrt
        figureOfMerit(iterIndex,1) = rms;
        
        % For triple cell
        % Calculate cc
%         figureOfMerit(iterIndex,2) = getPhaseCorrelation(FPhaseCorneredModel,FPhaseCorneredEstimate,FModulusCorneredModel,FModulusCorneredEstimate,HKLmask);
%       % Calculate R_factor
%         figureOfMerit(iterIndex,3) = getRFactor(FModulusCorneredModel,FModulusCorneredEstimate,HKLmask);
        
        % for unit cell
        % Calculate cc
        figureOfMerit(iterIndex,2) = getPhaseCorrelation(model_s_angle,objIter_s_phase,model_s_abs,objIter_s_abs,HKLmask);        
        % Calculate R_factor
        figureOfMerit(iterIndex,3) = getRFactor(model_s_abs,objIter_s_abs,HKLmask);
        % Calculate average phase_error
        figureOfMerit(iterIndex,4) = getaPhaseError(model_s_angle,objIter_s_phase,HKLmask);
        % Calculate weighted phase_error
        figureOfMerit(iterIndex,5) = getwPhaseError(model_s_angle,objIter_s_phase,model_s_abs,objIter_s_abs,HKLmask);

        if symConstraint == 1
            objEstimate = symC2221(objEstimate);
        end
        
        if mod(iterIndex,20) == 0
            fprintf('%dth iteration. rms: %f; CC: %f; R factor: %f \n',iterIndex,figureOfMerit(iterIndex,1),figureOfMerit(iterIndex,2),figureOfMerit(iterIndex,3));
        end
    end
    
    
    
    % Do Error reduction 10 times.
    for t = 1:iterER
        iterIndex = (i-1)*(iterHIO+iterER)+iterHIO + t+1;
        FComplexCorneredEstimatePrime = fftn(objEstimate);
%         FModulusCorneredEstimate = abs(FComplexCorneredEstimatePrime);
        
        % Current estimate of model in unit cell, |F| and phase
%         objIter = objEstimate(:,:,int32(sizeC/3+1):int32(sizeC*2/3));
        objIter = objEstimate;
        objIter_s_abs = abs(fftn(objIter));
        objIter_s_phase = angle(fftn(objIter));
        
        
        % The modulus constraint is implemented in the following line
        FComplexCorneredEstimated = FModulusCorneredModel.*exp(1i*(angle(FComplexCorneredEstimatePrime)));
        
        % Inverse Fourier Transform gives object.
        objEstimatePrime = ifftn(FComplexCorneredEstimated);
        
        % Impose support from real space
        objEstimate = abs(objEstimatePrime.*support1);
        
        % Calculate rms
        obj_outside_support = objEstimatePrime.*support2;
        rms = sqrt(sum(sum(sum(obj_outside_support.^2)))/sum(sum(sum(abs(objEstimatePrime.^2)))));
        figureOfMerit(iterIndex,1) = rms;

        % for triple cell
%         % Calculate cc
%         figureOfMerit(iterIndex,2) = getPhaseCorrelation(FPhaseCorneredModel,FPhaseCorneredEstimate,FModulusCorneredModel,FModulusCorneredEstimate,HKLmask);
%         
%         % Calculate R_factor
%         figureOfMerit(iterIndex,3) = getRFactor(FModulusCorneredModel,FModulusCorneredEstimate,HKLmask);
        
        % for unit cell
        % Calculate cc
        figureOfMerit(iterIndex,2) = getPhaseCorrelation(model_s_angle,objIter_s_phase,model_s_abs,objIter_s_abs,HKLmask);        
        % Calculate R_factor
        figureOfMerit(iterIndex,3) = getRFactor(model_s_abs,objIter_s_abs,HKLmask);
        % Calculate average phase_error
        figureOfMerit(iterIndex,4) = getaPhaseError(model_s_angle,objIter_s_phase,HKLmask);
        % Calculate weighted phase_error
        figureOfMerit(iterIndex,5) = getwPhaseError(model_s_angle,objIter_s_phase,model_s_abs,objIter_s_abs,HKLmask);
        
        if mod(iterIndex,20) == 0
            fprintf('%dth iteration. rms: %f; CC: %f; R factor: %f \n',iterIndex,figureOfMerit(iterIndex,1),figureOfMerit(iterIndex,2),figureOfMerit(iterIndex,3));
        end
        if symConstraint == 1
            objEstimate = symC2221(objEstimate);
        end
    end
    
    %  Perform shrinkwrap.
    if shrinkwrap == 1
        if guassianSigma > guassianSigmaMin
            guassianSigma = guassianSigma*(101-i)/100;
        else
            guassianSigma = guassianSigmaMin;
        end
        fprintf('Guassian sigma is %f: \n',guassianSigma);
        support1 = shrinkWrap(objEstimate,guassianSigma,guassianMaskThreshold);
        
        %  Should we impose hole constraint all the time?
        if constantSupportConstraint == 1
            support_1 = support_1.*support_0;
        end
        
        % Complimentary support.
        support2 = 1- support1; % complimentary support
    end
    
end

% obj = zeros([s_a,s_b,s_c]/3);

% Here support size is triple as the original one in each dimension.
% We just take the middle part, which is reconstructed object.
% objFinal = objEstimate(:,:,int32(sizeC/3+1):int32(sizeC*2/3)); 
objFinal = objEstimate;

objName = strcat('runID_',int2str(runID),'_obj.mat');
metricName = strcat('runID_',int2str(runID),'_metric.mat');
save(objName,'objFinal');
save(metricName,'figureOfMerit');

% Plot Figures
% RPhaseCCPlot(figureOfMerit,runID);
RPhaseCCPlot3(figureOfMerit,runID);

% Show reconstruction
show3DObjv2(objFinal,0.85,1,1);

% Show model
show3DObjv2(model,0.85,1,1);
