% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu 

function realLensCenters = FLFM_computeLensCenters(CaibrationImage, Camera, LensletGridModel)

% use the calibration image to retreive the transformations between the EI
transformationsStack = FLFM_retrieveEItransformations(LensletGridModel, CaibrationImage);

% lenslet centers according to the uniform spacing of the LensletGridModel
centersUniform = LFBuildGrid(LensletGridModel, Camera.gridType);

% actual lenslet centers (as from the EIs registration; the micro-lens are not perfectly uniform spaced in reality)
realLensCenters = zeros(size(centersUniform));
for j = 1 : size(transformationsStack,1)
    for k = 1 : size(transformationsStack,2)
        realLensCenters(j,k,1) = centersUniform(j,k,1) - transformationsStack{j,k}.T(3,1);
        realLensCenters(j,k,2) = centersUniform(j,k,2) - transformationsStack{j,k}.T(3,2);
%         fprintf('%.2f ',transformationsStack{j,k}.T(3,1:2));
    end
end

% make sure some centers are not outside the image (bad registration of the incomplete elemental images)
realLensCenters(abs(round(centersUniform) - round(realLensCenters)) > 50) = centersUniform(abs(round(centersUniform) - round(realLensCenters)) > 50);