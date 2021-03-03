% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu 

function Resolution = FLFM_computeResolution(LensletGridModel, Camera, depthRange, depthStep)

%% Compute sensor resolution
NspacingLenslet = [LensletGridModel.VSpacing, LensletGridModel.HSpacing]; % Number of pixels behind a lenslet 

% Corresponding sensor resolution
if(strcmp(Camera.gridType, 'hex'))
    sensorRes = [Camera.lensPitch*cosd(30)/NspacingLenslet(1), Camera.lensPitch./NspacingLenslet(2)];
    Nnum = [max(NspacingLenslet), max(NspacingLenslet) ];
end
if(strcmp(Camera.gridType, 'reg'))
    sensorRes = [Camera.lensPitch/NspacingLenslet(1), Camera.lensPitch./NspacingLenslet(2)];
    Nnum = NspacingLenslet;
end
Nnum = Nnum - (1-mod(Nnum,2));

%% Object space resolution (voxel size in um) 
texRes = sensorRes./Camera.M;
texRes(3) = depthStep;
% field of view in voxels
Resolution.fovRadVox = [round(Camera.fovRad./texRes(1)), round(Camera.fovRad./texRes(2))];
%% Set up a struct containing the resolution related info
Resolution.Nnum = Nnum;
Resolution.sensorRes = sensorRes;
Resolution.texRes = texRes;
Resolution.depthStep = depthStep;
Resolution.depthRange = depthRange;
Resolution.depths = depthRange(1) : depthStep : depthRange(2);
Resolution.superResFactor = Camera.superResFactor;