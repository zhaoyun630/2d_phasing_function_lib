function newObj = changeSolventFraction(obj,solventPercentage,solventDensity)

%**********************************************************************
% modificate solvent fraction by setting lower charges to constant values
%************************************************************************

[sortDensity, sortIndex] = sort(obj(:));
nVoxels = length(sortIndex);
indexPerc = floor(nVoxels*solventPercentage);
densityCutoff = sortDensity(indexPerc);


densityFilter = obj > densityCutoff;
newObj = obj.*densityFilter + (1-densityFilter)*solventDensity;
