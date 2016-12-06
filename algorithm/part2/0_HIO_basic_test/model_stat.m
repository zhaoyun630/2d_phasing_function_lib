% This piece of code is used to generate some statistic plot of model.

% read from file
% model = importdata('3rdu.mat');

% read from variable in workspace.
model = newObj;

% model label
modelLabel = '3rdx_1A';

model_ave = mean(model(:));
model_std = std(model(:));



% set number of bins
binNum = 100;

[counts,binValue] = hist(model(:),binNum);

% plot histogram 
figure(1)
plt1 = plot(binValue,counts);

% % save plot
% plotFile1 = strcat(modelLabel,'_hist.tif');
% saveas(plt1,plotFile1);

% plot cumCounts
cumCounts = cumsum(counts); % absolute cummulative counts
percCumCounts = cumCounts/cumCounts(binNum); % fractional cummulative counts
figure(2)
plot(binValue,percCumCounts);


% % save plot
% plotFile2 = strcat(modelLabel,'_cumDist.tif');
% saveas(plt1,plotFile2);