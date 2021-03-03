% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function [LensletGridModel, gridCoords] = LFM_processWhiteImage(WhiteImage, spacingPx, gridType, DebugBuildGridModel)

    GridModelOptions.ApproxLensletSpacing = spacingPx; %lensPitch / pixelPitch;
    GridModelOptions.Orientation = 'horz';
    GridModelOptions.FilterDiskRadiusMult = 1/3;
    GridModelOptions.CropAmt = 30; 
    GridModelOptions.SkipStep = 10;
    GridModelOptions.Precision = 'single';
    

    %---Find grid params---
    [LensletGridModel, gridCoords] = LFM_BuildLensletGridModel( WhiteImage, gridType, GridModelOptions, DebugBuildGridModel );