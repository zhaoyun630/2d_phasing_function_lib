function projectionPlot(obj,dimension,nDuplicate)

%*************************************************************************
% This function is used to calculate density projection on certain axis
% obj is 3D object. dimension can be set as 1, 2, 3, which means x,y,z.
% Editted by Yun Zhao on May 5th, 2015
%*************************************************************************

switch dimension
    case 1
        objZ = squeeze(sum(sum(obj,2),3));
        if nDuplicate >= 2
            objZ = repmat(objZ,[nDuplicate,1]);
        end
        layerNumber = linspace(1,length(objZ),length(objZ));
        plot(layerNumber,objZ);
    case 2
        objZ = squeeze(sum(sum(obj,1),3));
        if nDuplicate >= 2
            objZ = repmat(objZ,[nDuplicate,2]);
        end
        layerNumber = linspace(1,length(objZ),length(objZ));
        plot(layerNumber,objZ);
    case 3
        objZ = squeeze(sum(sum(obj,1),2));
        if nDuplicate >= 2
            objZ = repmat(objZ,[nDuplicate,3]);
        end
        layerNumber = linspace(1,length(objZ),length(objZ));
        plot(layerNumber,objZ);
end


