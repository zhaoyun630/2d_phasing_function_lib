% Used to obtain pointer support and object after solvent flattening.
% 
% Three types of pointer support:
% 1) Tight pointer support
% 2) Loose support from low threshold of experimental model
% 3) Tight pointer + random pointer
% 
% This code is used for testing the 1st case. The aim is to study what is
% the maximum a solvent can a cell have.
%
% The code is very long. 99% of them are useless. I am just too lazy to
% delete them. At the other point, I may use them just in case. There is no
% harm to leave it there.
% -Yun Zhao 2016.3.21

% Input parameter
obj = importdata('3rdx_c2221_2A.mat');
filterFrac = 0.1;
objThreshold = 0.6;
envThreshold = 0.36;
inputID = 6;

objSize = size(obj);
totalPixel = objSize(1)*objSize(2)*objSize(3);



%% Guassian filter
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
% GuassianFilter(:,:,1) = [0 0 0;0 1 0;0 0 0];
% GuassianFilter(:,:,2) = [0 1 0;1 1 1;0 1 0];
% GuassianFilter(:,:,2) = [0 1 0;1 20 1;0 1 0];
% GuassianFilter(:,:,3) = [0 0 0;0 1 0;0 0 0];
% GuassianFilter = GuassianFilter/13;

%% Set new object and create support

numSolventObj = round(objThreshold*totalPixel);
numSolventEnv = round(envThreshold*totalPixel);

sortDensity = sort(obj(:));
thresholdDensity = sortDensity(numSolventObj);
support1 = obj > thresholdDensity;
newObj = support1.*obj;
support2 = 1 - support1;

support3 = obj > sortDensity(numSolventEnv);

support4 = support3;

%% I don't understand why I wrote following code...
% randObj = rand(objSize).*support2;
% sortRandObj = sort(randObj(:));
% support3 = randObj > sortRandObj(totalPixel - numSolventEnv);


% % support3 = obj > thresholdEnvelope;
% 
% blurObj = newObj;
% 
% % support3 = support1;
% % i= 0
% % x1 = sum(support3(:))/totalPixel
% 
% blurObj = convn(blurObj,GuassianFilter,'same');
% sortEnvDensity = sort(blurObj(:));
% support4 = blurObj > sortEnvDensity(round(envThreshold*totalPixel));

% while x1<1-envThreshold
%     i = i+1
%     blurObj = convn(blurObj,GuassianFilter,'same');
% %     support3 = blurObj > minDensity;
% %     sortEnvDensity = sort(blurObj(:));
% %     support3 = blurObj > sortEnvDensity(round(envThreshold*totalPixel));
%      support3 = blurObj > 0;
%     x1 = sum(support3(:))/totalPixel
% end

envelope = support4;

% Generate triple support and obj.
[s1,s2,s3] = size(obj);
sizeA = 3*s1 -2;
sizeB = 3*s2 -2;
sizeC = 3*s3 -2;
triObj = zeros(s1,s2,sizeC);
triSupport = zeros(s1,s2,sizeC);
triObj(:,:,int32(sizeC/3+1):int32(sizeC*2/3))= newObj;
% triObj(:,:,int32(sizeC/3+1):int32(sizeC*2/3))= obj;
triSupport(:,:,int32(sizeC/3+1):int32(sizeC*2/3))= support4;

fcomplex = fftn(triObj);

modelFile = strcat('inputID',int2str(inputID),'_newObj.mat');
reflectionlist = strcat('inputID',int2str(inputID),'_fcomplex.mat');
supportFile = strcat('inputID',int2str(inputID),'_triSupport.mat');
triCellFile = strcat('inputID',int2str(inputID),'_triCell.mat');
envelopObjFile = strcat('inputID',int2str(inputID),'_envelopObj.mat');


% save(modelFile,'newObj');
% save(reflectionlist,'fcomplex');
% save(supportFile,'triSupport');
% % save(envelopObjFile,'blurObj');
% save(triCellFile,'triObj');

% calculate solvent fraction
solvent = 1 - sum(support4(:))/totalPixel

% calculate false constraint
a = support1 - support4;
b = (a>0);
support4 = 1-support4;
support_error = sum(b(:))/sum(support4(:))
