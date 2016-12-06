function [triSupport, triCell] = getTriCellSupport(newObj,envelope)

%% About this function
% This function pad zeros above and below to generate a triple cell and
% support

[sizeA, sizeB, sizeC] = size(newObj);
%  % Triple in c direction
triCell = zeros(sizeA,sizeB,3*sizeC-2);
triSupport = zeros(sizeA,sizeB,3*sizeC-2);
triCell(:,:,sizeC:2*sizeC-1) = newObj;
triSupport(:,:,sizeC:2*sizeC-1) = envelope;
% save('obj_abs_triple_3A.mat','triCell');
% save('support_triple_3A.mat','triSupport');





