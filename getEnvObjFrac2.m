function [envelope,newObj] = getEnvObjFrac2(obj,filterFrac,envThreshold,objThreshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a envelope support from object
% Input obj is a 3D array
% filterFrac is the fraction of half-length of unit cell in pixels. It
% ranges from 0 to 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Difference between getEnvObjFrac.m and getEnvObjFrac2.m

% getEnvObjFrac.m procedure:
% 1, Set a threshold for original object and find a support. Then
% use this support and original object to get a new object.
% 2, Set a lower threshold for original object and find a support.
% 3, Then apply guassian filter to this support. Then set an another
% threshold to get a final support for our algorithm.

% getEnvObjFrac2.m procedure:
% 1, Use guassian filter to blur the original object.
% 2, Set threshold for blurred object to get a support.
% 3, Set a higher threshold for blurred object to get an another
% support. And then use this support to filter original object to get a new
% obj.



%%

%% Part 1 Code for 3D Guassian filter

objSize = size(obj);

% because in reciprocal space, the size in each dimension is always
% odd(0,+/-1;+/-2;...) and symmetric. So we need to substract 1 and divide
% by 2 to get radius at certain resolution ring.
% pixelRes = resMax/round(((numPixel-1)/2));

numPixel = (objSize -1)/2;
filterSize = round(filterFrac*numPixel); % half of the filter size.

filterCoordX = -filterSize(1):1:filterSize(1);
filterCoordY = -filterSize(2):1:filterSize(2);
filterCoordZ = -filterSize(3):1:filterSize(3);

[x,y,z] = meshgrid(filterCoordX,filterCoordY,filterCoordZ);
GuassianFilterPrime = normpdf(x,0,filterSize(1)).*normpdf(y,0,filterSize(2)).*normpdf(z,0,filterSize(3));
GuassianFilter = GuassianFilterPrime/sum(sum(sum(GuassianFilterPrime)));

% % Now convolute the object density with a guassian filter
objConvGuass = imfilter(obj,GuassianFilter);

%% Part 2 Get envelope for unit cell

sortDensity = sort(objConvGuass(:));
thresholdDensity = sortDensity(round(envThreshold*objSize(1)*objSize(2)*objSize(3)));
envelope = objConvGuass > thresholdDensity;
% envelopeFile = strcat('envFrac',int2str(int32(filterFrac*100)),'picometer','.mat');
% save(envelopeFile,'envelope');

%% Part 3 Get modified object for unit cell
support2 = objConvGuass > sortDensity(round(objThreshold*objSize(1)*objSize(2)*objSize(3)));
newObj = support2.*obj;
% newObjFile = strcat('unitcellSig',int2str(int32(filterFrac*100)),'picometer','.mat');
% save(newObjFile,'newObj');
