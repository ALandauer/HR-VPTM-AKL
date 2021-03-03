% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function [H, Ht] = LFM_computeLFMatrixOperators(Camera, Resolution, LensletCenters)

%% Maximum spread area on the sensor and corresponding lens centers
% psf area is maximal at max depth
% compute psf size in lenslets, and image size in px
[~, p] = max(abs(Resolution.depthRange + Camera.offsetFobj)); % the PSF area is highest at furthest depth
maxDepth = Resolution.depthRange(p) + Camera.offsetFobj; %% for 2.0 add offest from f objective; center depths around Camera.dof
PSFsize = LFM_computePSFsize(maxDepth, Camera); 

% extract only those lenslets touched by the psf extent
usedLensletCenters = LFM_getUsedCenters(PSFsize, LensletCenters);
Resolution.usedLensletCenters = usedLensletCenters;

%% Object and ML space

IMGsizeHalf = max( Resolution.Nnum(2)*(PSFsize), 2*Resolution.Nnum(2)); % PSF size in pixels
disp(['Size of PSF IMAGE = ' num2str(IMGsizeHalf*2+1) 'X' num2str(IMGsizeHalf*2+1) ' [Pixels]']);

% yspace/xspace represent sensor coordinates
% yMLspace/xMLspace represent local ulens coordinates
Resolution.yspace = Resolution.sensorRes(1)*[-IMGsizeHalf:1:IMGsizeHalf];
Resolution.xspace = Resolution.sensorRes(2)*[-IMGsizeHalf:1:IMGsizeHalf];
Resolution.yMLspace = Resolution.sensorRes(1)* [- Resolution.Nnum_half(1) + 1 : 1 : Resolution.Nnum_half(1) - 1];
Resolution.xMLspace = Resolution.sensorRes(2)* [- Resolution.Nnum_half(2) + 1 : 1 : Resolution.Nnum_half(2) - 1];

%% Compute Patterns

% setup parallel pool
pool = gcp('nocreate');
if isempty(pool)
    parpool;
end

% compute native plane PSF for every depth
psfWaveStack = LFM_calcPSFAllDepths(Camera, Resolution);

tolLFpsf = 0.005; % clap small values in forward patterns to speed up computations
if (strcmp(Camera.focus,'single'))
    [H, Ht] = LFM_computePatternsSingleWaves(psfWaveStack, Camera, Resolution, tolLFpsf);
end

if (strcmp(Camera.focus,'multi'))
    [H, Ht] = LFM_computePatternsMultiWaves(psfWaveStack, Camera, Resolution, tolLFpsf);
end