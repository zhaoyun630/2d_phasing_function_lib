function originShift = getOrigin(obj,model)

%*************************************************************************
% This function is used to calculate origin shift by comparing model and
% structure which is reconstructed from iterative algorithms.
% Note: model and obj should share the same size.
% Created by Yun Zhao on May 18th, 2015
% Last editted by Yun Zhao on May 19th, 2015
%
% Comment: It's great, efficient!
%          It can be modified to more delicated algorithm, for example, 
%          searching couple of alternative vectors if there are many
%          symmetry operation elements.
%          Done!
%*************************************************************************

% Calculate projection to one axis for obj.
objX = squeeze(sum(sum(obj,2),3));
objY = squeeze(sum(sum(obj,1),3));
objZ = squeeze(sum(sum(obj,1),2));

% Calculate projection to one axis for model.
modelX = squeeze(sum(sum(model,2),3));
modelY = squeeze(sum(sum(model,1),3));
modelZ = squeeze(sum(sum(model,1),2));

% Calculate correlation.
[sizeX,sizeY,sizeZ] = size(obj);
corrX = zeros(sizeX,1);
corrY = zeros(sizeY,1);
corrZ = zeros(sizeZ,1);

for i=1:sizeX
    objXshift = circshift(objX,[i,0]);
    corrX(i) = dot(objXshift,modelX);
end

for j=1:sizeY
    objYshift = circshift(objY,[0,j]); 
    % Here I use [0 j] rather than [j 0] because it is a row, rather then column.
    corrY(j) = dot(objYshift,modelY);
end

for k=1:sizeZ
    objZshift = circshift(objZ,[k,0]);
    corrZ(k) = dot(objZshift,modelZ);
%     check where is the bug
%     fprintf('k is %d; and corrZ(k) is %d \n',k,corrZ(k));
%     fprintf('The first 5 number of objZshift is %d %d %d %d %d \n',...
%         objZshift(1),objZshift(2),objZshift(3),objZshift(4),objZshift(5));
%     pause
end


[~,xMaxPosition]=max(corrX);
[~,yMaxPosition]=max(corrY);
[~,zMaxPosition]=max(corrZ);
originShift = [xMaxPosition, yMaxPosition, zMaxPosition];



% figures are used to check Bugs
x=linspace(1,sizeX,sizeX);
y=linspace(1,sizeY,sizeY);
z=linspace(1,sizeZ,sizeZ);
figure(1)
subplot(2,1,1);
plot(x,objX,'b-*',x,modelX,'r-o')
subplot(2,1,2);
plot(x,corrX,'b-x')

figure(2)
subplot(2,1,1);
plot(y,objY,'b-*',y,modelY,'r-o')
subplot(2,1,2);
plot(y,corrY,'b-x')

figure(3)
subplot(2,1,1);
plot(z,objZ,'b-*',z,modelZ,'r-o')
subplot(2,1,2);
plot(z,corrZ,'b-x')


