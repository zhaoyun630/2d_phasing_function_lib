function show3DObjv1(obj,ShowPercent)

%*****************************************************************
% Used to display density reconstruction from Charge flipping.
% obj: any 2D or 3D array
% ShowPercent: range from 0~1. Within this fraction, voxels are not shown.
%               recommand starting with 0.97 for proteins.
% nReplicate: 0 or 1.replicate object with 2 copies in each dimension or not.
% degradeLevel: range 0~1. Cut the resolution.
%*****************************************************************


[objDensitySort,objDensitySortIndex] = sort(obj(:));
I_contour = objDensitySort(floor(ShowPercent*length(objDensitySortIndex)));

% Show at a fraction of original resolution.
% if degradeLevel == 0
%     objShow = obj;
% else
%     gObj = fftshift(fftn(obj));
%     gObjCenter = floor((size(gObj)+1)/2);
%     gObjHalfSize = floor((size(gObj)-1)/2);
%     gObjShowHalfSize = floor(gObjHalfSize*degradeLevel);
%     ShowIndexLow = int32(gObjCenter-gObjShowHalfSize);
%     ShowIndexHigh = int32(gObjCenter+gObjShowHalfSize);
%     gObjShow = gObj(ShowIndexLow(1):ShowIndexHigh(1),ShowIndexLow(2):ShowIndexHigh(2),ShowIndexLow(3):ShowIndexHigh(3));
%     objShow = ifftn(ifftshift(gObjShow));
% end
% 
% if nReplicate <= 1
%     objShow = objShow;
% else
%     objShow = repmat(objShow,[nReplicate 1 1]);
%     objShow = repmat(objShow,[1 nReplicate 1]);
%     objShow = repmat(objShow,[1 1 nReplicate]);
% end
objShow = obj;

    
[DimASize,DimBSize,DimCSize]=size(objShow);

figure()
fv = isosurface(objShow,I_contour);
p=patch(fv);
set(p,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
view(3); % axis tight
axis([1 DimASize 1 DimBSize 1 DimCSize])
% axis equal
camlight 
lighting gouraud
box on
xlabel('a');
ylabel('b');
zlabel('c');
