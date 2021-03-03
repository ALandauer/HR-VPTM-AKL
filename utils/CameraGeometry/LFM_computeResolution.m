% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function Resolution = LFM_computeResolution(LensletGridModel, TextureGridModel, Camera, depthRange, depthStep)

%% Compute sensor resolution
% Number of pixels behind a lenslet 
NspacingLenslet = [LensletGridModel.VSpacing, LensletGridModel.HSpacing];
NspacingTexture = [TextureGridModel.VSpacing, TextureGridModel.HSpacing];

% Corresponding sensor/tex resolution
if(strcmp(Camera.gridType, 'hex'))
    sensorRes = [Camera.lensPitch*cosd(30)/NspacingLenslet(1), Camera.lensPitch./NspacingLenslet(2)];
    Nnum = [max(NspacingLenslet) + 1, max(NspacingLenslet) + 1];
    TexNnum = [max(NspacingTexture) + 1, max(NspacingTexture) + 1];
end

if(strcmp(Camera.gridType, 'reg'))
    sensorRes = [Camera.lensPitch/NspacingLenslet(1), Camera.lensPitch./NspacingLenslet(2)];
    Nnum = NspacingLenslet;
    TexNnum = NspacingTexture;
end

Nnum = Nnum + (1-mod(Nnum,2));
TexNnum = TexNnum + (1-mod(TexNnum,2));

%% Size of a voxel in micrometers. (superResFactor is a factor of lensletResolution)
% When superResFactor == 1, we reconstruct at lenslet resolution
% When superResFactor == Nnum, we reconstruct at sensor resolution
% texScaleFactor = superResFactor./Nnum;

% make sure the superResFactor produces an odd number of voxels per
% repetition patch (in front of a mlens)
% TexNnum = floor(texScaleFactor(1)*Nnum);
% TexNnum = TexNnum + (1-mod(TexNnum,2));
texScaleFactor = TexNnum./Nnum;

texRes = sensorRes./(texScaleFactor*Camera.M);
texRes(3) = depthStep;


%% Compute mask for sensor/texture to avoid overlap during convolution (for hex grid)
% mask for the patches behind a lenslet 
sensMask = LFM_computePatchMask(NspacingLenslet, Camera.gridType, sensorRes, Camera.uRad, Nnum);

% texture mask (different of sensor mask when tex/image resolutions are decoupled)
texMask = LFM_computePatchMask(NspacingTexture, Camera.gridType, texRes, TexNnum(1)*texRes(1)/2, TexNnum);

%% Set up a struct containing the resolution related info
Resolution.Nnum = Nnum;
Resolution.Nnum_half = ceil(Nnum/2);
Resolution.TexNnum = TexNnum;
Resolution.TexNnum_half = ceil(TexNnum/2);
Resolution.sensorRes = sensorRes;
Resolution.texRes = texRes;
Resolution.sensMask = sensMask;
Resolution.texMask = texMask;
Resolution.depthStep = depthStep;
Resolution.depthRange = depthRange;
Resolution.depths = depthRange(1) : depthStep : depthRange(2);
Resolution.texScaleFactor = texScaleFactor;
Resolution.maskFlag = Camera.uLensMask;
Resolution.NspacingLenslet = NspacingLenslet;
Resolution.NspacingTexture = NspacingTexture;