% Resize a 3-dimensional object.

% Show at a fraction of original resolution.
degradeLevel = 0.7;
obj = importdata('inputID1_obj.mat');

if degradeLevel == 0
    objShow = obj;
else
    gObj = fftshift(fftn(obj));
    gObjCenter = floor((size(gObj)+1)/2);
    gObjHalfSize = floor((size(gObj)-1)/2);
    gObjShowHalfSize = floor(gObjHalfSize*degradeLevel);
    ShowIndexLow = int32(gObjCenter-gObjShowHalfSize);
    ShowIndexHigh = int32(gObjCenter+gObjShowHalfSize);
    gObjShow = gObj(ShowIndexLow(1):ShowIndexHigh(1),ShowIndexLow(2):ShowIndexHigh(2),ShowIndexLow(3):ShowIndexHigh(3));
    objShow = ifftn(ifftshift(gObjShow));
end