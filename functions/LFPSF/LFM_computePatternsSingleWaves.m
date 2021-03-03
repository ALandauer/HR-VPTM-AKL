% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function [H, Ht] = LFM_computePatternsSingleWaves(psfWaveStack, Camera, Resolution, tolLFpsf)

%% Precompute the MLA transmittance funtion
ulensPattern = LFM_ulensTransmittance(Camera, Resolution);
MLARRAY = LFM_mlaTransmittance(Camera, Resolution, ulensPattern);

%% Compute forward light trasport patterns -> patters for all texture points
disp('Computing forward patterns for single focus array')
H = LFM_computeForwardPatternsWaves(psfWaveStack, MLARRAY, Camera, Resolution);
H = ignoreSmallVals(H, tolLFpsf);

%% Compute backward light transport patterns -> patters for all sensor points
disp('Computing backward patterns for single focus array')
Ht = LFM_computeBackwardPatterns(H, Resolution, Camera.range, []); 
Ht = normalizeHt(Ht);
return