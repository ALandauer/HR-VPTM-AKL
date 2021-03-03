% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function lensCurrentType = LFM_extractLensType(lensletCenters, type)

lensVoxY = lensletCenters.vox(:,:,1);
lensVoxX = lensletCenters.vox(:,:,2);
lensPxY = lensletCenters.px(:,:,1);
lensPxX = lensletCenters.px(:,:,2);

% extract only those centers corresponding to lenslets of type "type"
lensTypes = lensletCenters.px(:,:,3);
lensCurrentType.px(:,:,1) = lensPxY(lensTypes == type);
lensCurrentType.px(:,:,2) = lensPxX(lensTypes == type);
lensCurrentType.vox(:,:,1) = lensVoxY(lensTypes == type);
lensCurrentType.vox(:,:,2) = lensVoxX(lensTypes == type);
