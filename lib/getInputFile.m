function getInputFile()

%% Explanation about this function
% This function is supposed to take an input txt file as input.The input 
% file saves all your simulation conditions. But at this stage, we 
% just give all initial conditions inside the funciton. 

%% Part 1 Input parameters.
hklFile = '3rdx_full_3A.hkl';
inputID = 3;
densityShift = 1; % it means to shift the overall charges so that the minimum charge is 0.
filterFrac = 0.1; % 
envThreshold = 0.4; % lower value means less solvent fraction
objThreshold = 0.8; % lower value means less solvent fraction
envMethod = 1;

inputPara = strcat('inputID',int2str(inputID),'_inputPara.mat');
save(inputPara,'hklFile','inputID','filterFrac', ...
    'densityShift','envThreshold','objThreshold','envMethod');
%% Part 2 Running
%% Get object from hkl File.
fprintf('Calculating density map from HKL file ...\n');
obj = hkl2obj(hklFile, densityShift);
objFile = strcat('inputID',int2str(inputID),'_obj.mat');
save(objFile,'obj');

% Get an molecular envelope and a modified object.
fprintf('Calculating molecular envelope and modify object ... \n');
if envMethod == 1
    [envelope,newObj] = getEnvObjFrac(obj,filterFrac,envThreshold,objThreshold);
elseif envMethod == 2
    [envelope,newObj] = getEnvObjAbs(obj,filterFrac,envThreshold,objThreshold);
elseif envlope == 3
    % The following function hasn't been implemented yet
    % It is supposed to use getEnvObjStd().
    [envelope,newObj] = getEnvObjFrac(obj,filterFrac,envThreshold,objThreshold); 
else
    fprintf('Invalid method. Set envMethod to 1, 2 or 3.\n');
end


envelopeFile = strcat('inputID',int2str(inputID),'_envelope.mat');
newObjFile = strcat('inputID',int2str(inputID),'_newObj.mat');

save(envelopeFile,'envelope');
save(newObjFile,'newObj');


%% Pad envelope and newobj to a triple cell.
fprintf('Pad triple cell and support ... \n');
[triSupport,triCell] = getTriCellSupport(newObj,envelope);
triSupportFile = strcat('inputID',int2str(inputID),'_triSupport.mat');
triCellFile = strcat('inputID',int2str(inputID),'_triCell.mat');
save(triCellFile,'triCell');
save(triSupportFile,'triSupport');

%% Generate moludus from triple cell.
% Note here the F(000) is at the corner of 3D array.
fprintf('Generating complex structure factor from triple cell ... \n');

fcomplex = fftn(triCell);
fcomplexFile = strcat('inputID',int2str(inputID),'_fcomplex.mat');
save(fcomplexFile,'fcomplex');
% modulus = abs(fftn(triCell));
% phases = angle(fftn(triCell));
% modulusFile = strcat('inputID',int2str(inputID),'_modulus.mat');
% phaseFile = strcat('inputID',int2str(inputID),'_phases.mat');
% save(modulusFile,'modulus');
% save(phaseFile,'phases');