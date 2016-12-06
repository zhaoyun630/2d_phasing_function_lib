% plot density distribution

model = importdata('unit_cell.mat');
obj = importdata('gkxy.mat');
nDuplicate = 2;

% if obj is triple cell then use the following code
[sizeA,sizeB,sizeC] = size(obj);
obj1 = obj(:,:,int32(sizeC/3+1):int32(sizeC*2/3));
obj = obj1;

% plot along x
objZ1 = squeeze(sum(sum(obj,2),3));
model1 = squeeze(sum(sum(model,2),3));
if nDuplicate >= 2
    objZ1 = repmat(objZ1,[nDuplicate,1]);
    model1 = repmat(model1,[nDuplicate,1]);
end
layerNumber1 = linspace(1,length(objZ1),length(objZ1));
subplot(3,1,1);
plot(layerNumber1,objZ1,'r',layerNumber1,model1,'b');
xlabel('a','fontsize',15);
% ylabel('integrated charge density');
legend('HIO estimate','model');
title('integrated charge density distribution over two unit cells','fontsize',18)


% plot along y
objZ2 = squeeze(sum(sum(obj,1),3));
model2 = squeeze(sum(sum(model,1),3));
if nDuplicate >= 2
    objZ2 = repmat(objZ2,[1,nDuplicate]);
    model2 = repmat(model2,[1,nDuplicate]);
end
layerNumber2 = linspace(1,length(objZ2),length(objZ2));
subplot(3,1,2);
plot(layerNumber2,objZ2,'r',layerNumber2,model2,'b');
xlabel('b','fontsize',15);
ylabel('integrated charge density','fontsize',18);
legend('HIO estimate','model');

% plot along z
objZ3 = squeeze(sum(sum(obj,1),2));
model3 = squeeze(sum(sum(model,1),2));
if nDuplicate >= 2
    objZ3 = repmat(objZ3,[nDuplicate,1]);
    model3 = repmat(model3,[nDuplicate,1]);
end
layerNumber3 = linspace(1,length(objZ3),length(objZ3));
subplot(3,1,3);
plot(layerNumber3,objZ3,'r',layerNumber3,model3,'b');
xlabel('c','fontsize',15);
legend('HIO estimate','model');
% ylabel('integrated charge density');

