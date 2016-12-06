function [density1,accumulateProbability] = densityStat(obj,percentage,nBins)

%**********************************************************************
% Firstly, it may be used to plot denisty distribution
% Secondly, it calculate the accumulate probability at a given percent.
%************************************************************************

[sortDensity, sortIndex] = sort(obj(:));
nVoxels = length(sortIndex);
xVoxel = linspace(1,nVoxels,nVoxels);
indexPerc = floor(nVoxels*percentage);
density1 = sortDensity(indexPerc);
minDensity = min(sortDensity);
maxDensity = max(sortDensity);

densityRange = linspace(minDensity,maxDensity,(nBins+1));
densityCounts = zeros((nBins+1),1);
densityInterval = (maxDensity-minDensity)/nBins;

for i=1:nVoxels
    intervalIndex = floor((sortDensity(i)-minDensity)/densityInterval)+1;
    densityCounts(intervalIndex) = densityCounts(intervalIndex)+1;
end

% if plotfigure == 1
figure()
    subplot(2,1,1)
    
    plot(densityRange,densityCounts,'r');
    xlabel('charge density','FontSize', 15,'FontWeight','bold');
    ylabel('counts','FontSize', 15,'FontWeight','bold');
    title('charge density distribution within unit cell','FontSize', 15,'FontWeight','bold');
% elseif plotfigure ==2
    subplot(2,1,2)
    plot(xVoxel,sortDensity,'m');
    xlabel('pixel number','FontSize', 15,'FontWeight','bold');
    ylabel('charge density','FontSize', 15,'FontWeight','bold');
    title('sorted charge density over voxel number','FontSize', 15,'FontWeight','bold');
% else
% end

accumulateProbability = sum(densityCounts(1:floor((density1-minDensity)/densityInterval)))/nVoxels;