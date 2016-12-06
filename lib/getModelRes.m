% This script is used to generate a model at a lower resoluton with the
% same 3D array size as original model.

function obj = getModelRes(modelFile,degradeLevel)

% *****************************************************************
% modelFile: model calculated from high resolution hkl list
% degradeLevel: the fraction of hkl values which are used to calculate new
% model
% How does this piece of code work?
% It first perform FFT on original model; then certain fraction of low
% resolution hkl values are taken from the complete hkl list. Then zeros
% are padded in outter shell to make its size the same as original one. At
% last, the new model is calculated from new hkl list.

% Created by Yun Zhao on April 27,2016


model = importdata(modelFile);
hklArray = fftshift(fftn(model));
hklNew = zeros(size(hklArray));

hklCenter = floor((size(hklArray)+1)/2);
hklHalfSize = floor((size(hklArray)-1)/2);
hklShowHalfSize = floor(hklHalfSize*degradeLevel);
ShowIndexLow = int32(hklCenter-hklShowHalfSize);
ShowIndexHigh = int32(hklCenter+hklShowHalfSize);
hklShow = hklArray(ShowIndexLow(1):ShowIndexHigh(1),ShowIndexLow(2):ShowIndexHigh(2),ShowIndexLow(3):ShowIndexHigh(3));
hklNew(ShowIndexLow(1):ShowIndexHigh(1),ShowIndexLow(2):ShowIndexHigh(2),ShowIndexLow(3):ShowIndexHigh(3)) = hklShow;
obj = abs(ifftn(ifftshift(hklNew)));
