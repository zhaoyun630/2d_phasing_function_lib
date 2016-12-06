% Used to obtain support and object after solvent flattening.

% % Method:
% 	1) apply a Gaussian filter for high resolution pdb model to get an envelope S1 with certain solvent constant.
% 	2) apply the same Gaussian filter to low resolution pdb model to get another envelope S2 with same solvent constant.
% 	3) use the model inside S1 to generate structure factor amplitudes A1. The model is denoted as MS1.
% 	4) use S2 and A1 to reconstruct MS1.

% Attention: This approach used here is very different from getEnvObj.m in
% folder /simulation_thesis_part1/.



% Input parameter
obj_model = importdata('3rdx.mat');
obj_mr = importdata('3rdx_10.mat'); % model used for Molecular replacement.
filterFrac = 0.2;
objThreshold = 0.52;
envThreshold = 0.50;
inputID = 2;

objSize = size(obj_model);
totalPixel = objSize(1)*objSize(2)*objSize(3);

numPixel = (objSize -1)/2;
filterSize = round(filterFrac*numPixel); % half of the filter size.

for i=1:3
    if filterSize(i)<1
        filterSize(i)=1;
    end
end

guassianSigma = 1; % The unit is the half length of Guassian filter.

filterCoordX = -filterSize(1):1:filterSize(1);
filterCoordY = -filterSize(2):1:filterSize(2);
filterCoordZ = -filterSize(3):1:filterSize(3);



[x,y,z] = meshgrid(filterCoordX,filterCoordY,filterCoordZ);
GuassianFilterPrime = normpdf(x,0,guassianSigma*filterSize(1)).* ...
    normpdf(y,0,guassianSigma*filterSize(2)).*normpdf(z,0,guassianSigma*filterSize(3));
GuassianFilter = GuassianFilterPrime/sum(sum(sum(GuassianFilterPrime)));


blurObj = obj_model;
% blurMr = obj_mr;
blurMr = obj_model;

% support3 = support1;
% i= 0
% x1 = sum(support3(:))/totalPixel

blurObj = convn(blurObj,GuassianFilter,'same');
blurMr = convn(blurMr,GuassianFilter,'same');

sortEnvDensity = sort(blurObj(:));
support1 = blurObj > sortEnvDensity(round(objThreshold*totalPixel));

sortEnvDensity_Mr = sort(blurMr(:));
support2 = blurMr > sortEnvDensity_Mr(round(envThreshold*totalPixel));


% while x1<1-envThreshold
%     i = i+1
%     blurObj = convn(blurObj,GuassianFilter,'same');
% %     support3 = blurObj > minDensity;
% %     sortEnvDensity = sort(blurObj(:));
% %     support3 = blurObj > sortEnvDensity(round(envThreshold*totalPixel));
%      support3 = blurObj > 0;
%     x1 = sum(support3(:))/totalPixel
% end

newObj = support1.*obj_model;
envelope = support2;

% Generate triple support and obj.
[s1,s2,s3] = size(obj_model);
sizeA = 3*s1 -2;
sizeB = 3*s2 -2;
sizeC = 3*s3 -2;
triObj = zeros(s1,s2,sizeC);
triSupport = zeros(s1,s2,sizeC);
triObj(:,:,int32(sizeC/3+1):int32(sizeC*2/3))= newObj;
triSupport(:,:,int32(sizeC/3+1):int32(sizeC*2/3))= envelope;

fcomplex = fftn(triObj);

modelFile = strcat('inputID',int2str(inputID),'_newObj.mat');
reflectionlist = strcat('inputID',int2str(inputID),'_fcomplex.mat');
supportFile = strcat('inputID',int2str(inputID),'_triSupport.mat');
triCellFile = strcat('inputID',int2str(inputID),'_triCell.mat');
envelopObjFile = strcat('inputID',int2str(inputID),'_envelopObj.mat');


save(modelFile,'newObj');
save(reflectionlist,'fcomplex');
save(supportFile,'triSupport');
save(envelopObjFile,'blurMr');
save(triCellFile,'triObj');

% calculate solvent fraction
solvent = 1 - sum(support1(:))/totalPixel

% calculate false constraint
a = support2 - support1;
b = (a<0);
support3 = 1-support1;
c = sum(b(:))/sum(support3(:))
