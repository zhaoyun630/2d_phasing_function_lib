% ########################################################################
%            shrinkwrap function for 3D object
% ########################################################################
%
% First created by Yun Zhao on May 26th, 2015
%   
%*************************************************************************
% It takes Guassian filter parameters and obj as input. And it will return
% a new support for iterative phasing algorithm.

% Editted again by Yun Zhao on Sep 17th, 2015
% I replaced imfilter with convn. imfilter doesn't work in a way I
% expected. I should look into imfilter in detail some time.
%*************************************************************************

function newSupport = shrinkWrap(objEstimate,guassianSigma,guassianThreshold)


%  guassianFilterSize = guassianSigma;
maskCoord= -guassianSigma:1:guassianSigma;
[gmask_x,gmask_y,gmask_z] = meshgrid(maskCoord,maskCoord,maskCoord);
GuassianFilterPrime = normpdf(gmask_x,0,guassianSigma).*normpdf(gmask_y,0,guassianSigma).*normpdf(gmask_z,0,guassianSigma);
GuassianFilter = GuassianFilterPrime/sum(sum(sum(GuassianFilterPrime)));

% Now convolute the object density with a guassian filter
objEstimateConvGuass = convn(objEstimate,GuassianFilter,'same');

% Calculate some statistic numbers for support update
newSupport = objEstimateConvGuass > guassianThreshold*max(objEstimateConvGuass(:));

