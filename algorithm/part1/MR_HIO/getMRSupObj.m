% Used to obtain Molecular Replacement(MR) support and object after solvent flattening.
% 
% MR support description:
% 1) First get a tight pointer support
% 2) Set a box region to 1 where ligands are supposed to bind.
% -Yun Zhao 2016.3.22

% Input parameter
% Experimental Structure
obj = importdata('inputID1_obj.mat');

% Set solvent fraction
objThreshold = 0.45;

% Set block size
s1 = 10;
s2 = 10;
s3 = 16;
ligandBlock = ones(s1,s2,s3);

% Set block location 
% Fractional coordinate, which ranges from 0 to 1.
f1 = 0.14;
f2 = 0.14;
f3 = 0.3;

% Label input files
inputID = 1;
%% Get a tight pointer support and new object
objSize = size(obj);
totalPixel = objSize(1)*objSize(2)*objSize(3);
numSolventObj = round(objThreshold*totalPixel);
sortDensity = sort(obj(:));
support1 = obj > sortDensity(numSolventObj); % Tight Point Mask
support1 = support1.*1.0;

newObj = support1.*obj;

%% Locate block
c1 = floor(f1*objSize(1))+1;
c2 = floor(f2*objSize(2))+1;
c3 = floor(f3*objSize(3))+1;

if c1+s1-1>objSize(1)
    error('Block does not fit support. Set a smaller f1 value.\n');
end

if c2+s2-1>objSize(2)
    error('Block does not fit support. Set a smaller f2 value.\n');
end

if c3+s2-1>objSize(3)
    error('Block does not fit support. Set a smaller f3 value.\n');
end

support2 = support1;

support2(c1:c1+s1-1,c2:c2+s2-1,c3:c3+s3-1)=ligandBlock;

%% Create a demo image for block
demoObj = newObj;
demoObj(c1:c1+s1-1,c2:c2+s2-1,c3:c3+s3-1)=ligandBlock*max(newObj(:));

%% Generate triple support and obj.
[s1,s2,s3] = size(obj);
sizeA = 3*s1 -2;
sizeB = 3*s2 -2;
sizeC = 3*s3 -2;
triObj = zeros(s1,s2,sizeC);
triSupport = zeros(s1,s2,sizeC);
triObj(:,:,int32(sizeC/3+1):int32(sizeC*2/3))= newObj;
triSupport(:,:,int32(sizeC/3+1):int32(sizeC*2/3))= support2;

fcomplex = fftn(triObj);

%% Save input data
modelFile = strcat('inputID',int2str(inputID),'_newObj.mat');
reflectionlist = strcat('inputID',int2str(inputID),'_fcomplex.mat');
supportFile = strcat('inputID',int2str(inputID),'_triSupport.mat');
triCellFile = strcat('inputID',int2str(inputID),'_triCell.mat');
demoObjFile = strcat('inputID',int2str(inputID),'_demoObj.mat');


save(modelFile,'newObj');
save(reflectionlist,'fcomplex');
save(supportFile,'triSupport');
save(demoObjFile,'demoObj');
save(triCellFile,'triObj');

%% A double check of support
% calculate solvent fraction
solvent = 1 - sum(support2(:))/totalPixel

% calculate false constraint
a = support1 - support2;
b = (a>0);
support3 = 1-support2;
c = sum(b(:))/sum(support3(:))
