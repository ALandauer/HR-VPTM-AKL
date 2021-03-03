% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function [H, Ht] = LFM_computePatternsMultiWaves(psfWaveStack, Camera, Resolution, tolLFpsf)

%% Compute patters for every position in the ulens range. Hex grid discretization is not symmetric about the center of the lens.
disp('Computing forward patters for multifocal array')

tic
%% Pattern 1
ulensPattern = LFM_ulensTransmittance(Camera, Resolution);
MLARRAY = LFM_mlaTransmittance(Camera, Resolution, ulensPattern);
H1 = LFM_computeForwardPatternsWaves(psfWaveStack, MLARRAY, Camera, Resolution);

%% Pattern 2
CameraShift = Camera;
CameraShift.fm = circshift(Camera.fm, -1);
ulensPattern = LFM_ulensTransmittance(CameraShift, Resolution);
MLARRAY = LFM_mlaTransmittance(CameraShift, Resolution, ulensPattern);
H2 = LFM_computeForwardPatternsWaves(psfWaveStack, MLARRAY, CameraShift, Resolution);

%% Pattern 3
CameraShift.fm = circshift(Camera.fm, -2);
ulensPattern = LFM_ulensTransmittance(CameraShift, Resolution);
MLARRAY = LFM_mlaTransmittance(CameraShift, Resolution, ulensPattern);
H3 = LFM_computeForwardPatternsWaves(psfWaveStack, MLARRAY, CameraShift, Resolution);
toc   

%% ignore small values
H = cell(1, 3); 
H1 = ignoreSmallVals(H1, tolLFpsf);
H2 = ignoreSmallVals(H2, tolLFpsf);
H3 = ignoreSmallVals(H3, tolLFpsf);

H{1,1} = H1;
H{1,2} = H2;
H{1,3} = H3;

%% Compute backward light transport patterns -> patterns for all texture points
disp('Computing backward patters for multifocal array')
tic
lensOrder = [1,2,3];
Ht1 = LFM_computeBackwardPatterns(H, Resolution, Camera.range, lensOrder); 
Ht2 = LFM_computeBackwardPatterns(H, Resolution, Camera.range, circshift(lensOrder, -1)); 
Ht3 = LFM_computeBackwardPatterns(H, Resolution, Camera.range, circshift(lensOrder, -2)); 

Ht = cell(1, 3); 
Ht1 = normalizeHt(Ht1);
Ht2 = normalizeHt(Ht2);
Ht3 = normalizeHt(Ht3);

Ht{1,1}  = Ht1;
Ht{1,2}  = Ht2;
Ht{1,3}  = Ht3;
toc