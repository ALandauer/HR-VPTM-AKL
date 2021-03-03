% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu & Josue Page
%
% Updated for LFM2PT, AKL 2021-03-03
%

% This function imports the dependecies into the workspace
function import2ws()


addpath('functions/');

%Trial-MPT specific
addpath('utils/external/Scatter2Grid3D/');
addpath('utils/external/MPT_src/');

%oLaF specific
addpath('functions/projectionOperators/');
addpath('functions/LFPSF/');
addpath('functions/LFPSF/MLATransmittance/');
addpath('utils/CameraGeometry/');
addpath('utils/ImageRectification/');
addpath('utils/Aliasing/');
addpath('utils/misc/');
addpath('utils/misc/colormaps/');
addpath('utils/external/yaml_io/');


