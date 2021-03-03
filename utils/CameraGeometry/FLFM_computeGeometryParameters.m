% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu 

function [LensletCenters, Resolution] = FLFM_computeGeometryParameters(CaibrationImage, Camera, LensletGridModel, depthRange, depthStep)

%% Compute resolution according to the new grid
Resolution = FLFM_computeResolution(LensletGridModel, Camera, depthRange, depthStep);
disp(['Super resolution factor of: ', num2str(Resolution.superResFactor),' Pix size: [', num2str(Resolution.sensorRes(1)),', ',num2str(Resolution.sensorRes(2)),...
        '] Vox size: [', num2str(Resolution.texRes(1)), ', ',num2str(Resolution.texRes(2)),', ',num2str(Resolution.texRes(3)),']']);
    
%% Compute lenslets centers on the sensor and corresponding repetition patches centers in texture space (in voxels)
LensletCenters = FLFM_computeLensCenters(CaibrationImage, Camera, LensletGridModel); 
Resolution.LensletCenters = LensletCenters;

% sensor size in pixels
Resolution.sensorSize = size(CaibrationImage);