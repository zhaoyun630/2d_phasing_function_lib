function [envelope,newObj] = getEnvObjFrac3(obj,filterFrac,envThreshold,objThreshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a envelope support from object
% Input obj is a 3D array
% filterFrac is the fraction of half-length of unit cell in pixels. It
% ranges from 0 to 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Difference between getEnvObjFrac3.m and getEnvObjFrac2.m, get EnvObjFrac.m

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

% getEnvObjFrac3.m procedure:
% It will generate a correct support.
% 1, Set threashold for original object and find a tight support.
% 2, Apply Gaussian filter to the object and calculate a new filter
% 3, If the volume of new filter exceed expectation, stop; If not, iterate.


%%

%% Part 1 Code for 3D Guassian filter
% Resoluton is related to the size of pixels.
% The unit of sigma is pixels.
% resMax = 3; % Maximum resolution
objSize = size(obj);
totalPixel = objSize(1)*objSize(2)*objSize(3);

% because in reciprocal space, the size in each dimension is always
% odd(0,+/-1;+/-2;...) and symmetric. So we need to substract 1 and divide
% by 2 to get radius at certain resolution ring.
% pixelRes = resMax/round(((numPixel-1)/2));

numPixel = (objSize -1)/2;
filterSize = round(filterFrac*numPixel); % half of the filter size.

for i=1:3
    if filterSize(i)<1
        filterSize(i)=1;
    end
end

% fprintf('filterSize is: \n');
% filterSize
guassianSigma = 1; % The unit is the half length of Guassian filter.

filterCoordX = -filterSize(1):1:filterSize(1);
filterCoordY = -filterSize(2):1:filterSize(2);
filterCoordZ = -filterSize(3):1:filterSize(3);

[x,y,z] = meshgrid(filterCoordX,filterCoordY,filterCoordZ);
GuassianFilterPrime = normpdf(x,0,guassianSigma*filterSize(1)).* ...
    normpdf(y,0,guassianSigma*filterSize(2)).*normpdf(z,0,guassianSigma*filterSize(3));
GuassianFilter = GuassianFilterPrime/sum(sum(sum(GuassianFilterPrime)));

% % Now convolute the object density with a guassian filter
% objConvGuass = imfilter(obj,GuassianFilter);

%% Part 2 Get envelope for unit cell

sortDensity = sort(obj(:));
thresholdDensity = sortDensity(round(objThreshold*totalPixel));
support1 = obj > thresholdDensity;
newObj = support1.*obj;
% envelope = imfilter(support1,GuassianFilter) > 0;
% envelope = imfilter(support1,GuassianFilter);
support2 = support1;
x1 = sum(support2(:))/totalPixel;

i=0;
while x1>envThreshold
    newObj = convn(newObj,GuassianFilter,'same');
    support2 = newObj > 0;
    x1 = sum(support2(:))/totalPixel;
    i
end

envelope = support2;
% % Use multiple convolution on support. I though it would make support more
% % soomth, but it doesn't help.
% % for i=1:10
% %     envelope1 = convn(envelope1,GuassianFilter,'same');
% % end

% sortEnvelope1 = sort(envelope1(:));
% thresholdDensity1 = sortEnvelope1(round(envThreshold*totalPixel));
% envelope = envelope1 > thresholdDensity1;
% % envelope = envelope1 > 0;
% % envelope = support1;
% 
% % envelope = obj > thresholdDensity;
% 
% % envelopeFile = strcat('envFrac',int2str(int32(filterFrac*100)),'picometer','.mat');
% % save(envelopeFile,'envelope');

%% Part 3 Get modified object for unit cell
% support2 = obj > sortDensity(round(objThreshold*totalPixel));
% newObj = support2.*obj;
% % newObjFile = strcat('unitcellSig',int2str(int32(filterFrac*100)),'picometer','.mat');
% % save(newObjFile,'newObj');
