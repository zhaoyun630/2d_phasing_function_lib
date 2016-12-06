function [envelope,newObj] = getEnvObjAbs(obj,filterFrac,envThreshold,objThreshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a envelope support from object
% Input obj is a 3D array
% sigmaInRes is in unit with angstrom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%% Part 1 Code for 3D Guassian filter
% Resoluton is related to the size of pixels.
% The unit of sigma is pixels.
resMax = 3; % Maximum resolution
numPixel = max(size(obj));

% because in reciprocal space, the size in each dimension is always
% odd(0,+/-1;+/-2;...) and symmetric. So we need to substract 1 and divide
% by 2 to get radius at certain resolution ring.
pixelRes = resMax/round(((numPixel-1)/2));
sigmaInPixel = filterFrac/pixelRes;
filterSize = round(sigmaInPixel);
filterCoord = -filterSize:1:filterSize;

[x,y,z] = meshgrid(filterCoord,filterCoord,filterCoord);
GuassianFilterPrime = normpdf(x,0,sigmaInPixel).*normpdf(y,0,sigmaInPixel).*normpdf(z,0,sigmaInPixel);
GuassianFilter = GuassianFilterPrime/sum(sum(sum(GuassianFilterPrime)));

% % Now convolute the object density with a guassian filter
objConvGuass = imfilter(obj,GuassianFilter);

%% Part 2 Get envelope for unit cell
% % 1st way to get support: contour at = mean + std.
% densityMean = mean(objConvGuass(:));
% densitySTD = std(objConvGuass(:));
% envelope = objConvGuass > (densityMean+densitySTD);

% 2nd way to get support: contour at fraction*Max.
densityMax = max(objConvGuass(:));
envelope = objConvGuass > envThreshold*densityMax;
envelopeFile = strcat('envSig',int2str(int32(filterFrac*100)),'picometer','.mat');
% save(envelopeFile,'envelope');

% 3rd way is to cut at certain fraction of volume.

%% Part 3 Get modified object for unit cell
support2 = objConvGuass > objThreshold*densityMax;
newObj = support2.*obj;
newObjFile = strcat('unitcellSig',int2str(int32(filterFrac*100)),'picometer','.mat');
% save(newObjFile,'newObj');












%%
% % Part 2: Original code to do 3D Gaussian filter. Not recommended. Use code
% % in Part 1.
% % The following code difines filter size and Gaussian sigma separately
% % I think it is good for 3D protein case. We only want to manipulate
% % resolution, which gives sigma.
% % As to the filter size, there is no harm to choose a big one. 
% % So I commented out the code in below. Use the code in 1st part.

% maskSize = round(sigma*size(obj));
% guassianSigma = 1;
% maskCoordA = -maskSize(1):1:maskSize(1);
% maskCoordB = -maskSize(2):1:maskSize(2);
% maskCoordC = -maskSize(3):1:maskSize(3);
% [gmask_x,gmask_y,gmask_z] = meshgrid(maskCoordA,maskCoordB,maskCoordC);
% GuassianFilterPrime = normpdf(gmask_x,0,guassianSigma).*normpdf(gmask_y,0,guassianSigma).*normpdf(gmask_z,0,guassianSigma);
% GuassianFilter = GuassianFilterPrime/sum(sum(sum(GuassianFilterPrime)));
% 
% % Now convolute the object density with a guassian filter
% objConvGuass = imfilter(obj,GuassianFilter);
% support = objConvGuass > 0.1;

