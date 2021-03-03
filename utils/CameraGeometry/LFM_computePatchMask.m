% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function mask = LFM_computePatchMask(spacing, gridType, res, patchRad, Nnum)

ysensorspace = [(-floor(Nnum(1)/2)):1:floor(Nnum(1)/2)];
xsensorspace = [(-floor(Nnum(2)/2)):1:floor(Nnum(2)/2)];
[x,y] = meshgrid(res(1)*ysensorspace, res(2)*xsensorspace);
mask = sqrt(x.*x+y.*y) < patchRad;

% Resolve for holes and overlaps
mask = LFM_fixMask(mask, spacing, gridType);