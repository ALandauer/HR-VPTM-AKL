% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu 

function [Camera, LensletGridModel] = FLFM_setCameraParams(configFile, superResFactor)

%% set FLFM parameters:
%%%% microscope params
% fobj-> objective magnification
% NA-> objective aperture
% f1 -> relay lens 1
% f2 -> relay lens 2

%%%% sensor
% lensPitch-> lenslet pitch
% pixelPitch-> sensor pixel pitch

%%%% MLA params
% gridType-> microlens grid type: "reg" -> regular grid array; "hex" -> hexagonal grid array 
% fm-> focal length of the lenslets

%%%% light characteristics
% n-> refraction index (1 for air)
% wavelenght-> wavelenght of the the emission light

%%%% distances
% mla2sensor-> distance between MLA and sensor

%%% MLA array descriptor 
% spacingPixels-> number of pixels between horizontal neighboring elemental (sub-aperture) images
% noLensHoriz-> number of elemental images horizontally 
% noLensVert-> number of elemental images vertically
% shiftRow-> '1' or '2' when odd or even rows are shifted, respectively (in hex grid only)

Camera = ReadYaml(configFile);
Camera.spacingPixels = Camera.spacingPixels * superResFactor;
Camera.superResFactor = superResFactor;

Camera.objRad = Camera.fobj * Camera.NA; % objective radius
Camera.k = 2*pi*Camera.n/Camera.WaveLength; % wave number
Camera.M = Camera.fm*Camera.f1/(Camera.f2*Camera.fobj); %total system magnification

% field stop radius 
Camera.fsRad = Camera.lensPitch/2 * Camera.f2/Camera.fm;

% field of view radius
Camera.fovRad =  Camera.fsRad * Camera.fobj/Camera.f1;

%%%% MLA array descriptor
LensletGridModel.gridType = Camera.gridType;
LensletGridModel.UMax = Camera.noLensHoriz;
LensletGridModel.VMax = Camera.noLensVert;
LensletGridModel.FirstPosShiftRow = Camera.shiftRow; % in hexagonal grids 
LensletGridModel.Orientation = 'horz';

LensletGridModel.HSpacing = Camera.spacingPixels;
LensletGridModel.HSpacing = LensletGridModel.HSpacing + mod(LensletGridModel.HSpacing,2); % make sure it is even

if(strcmp(Camera.gridType, 'hex'))
    LensletGridModel.VSpacing = round(sqrt(3)/2*LensletGridModel.HSpacing);
    LensletGridModel.VSpacing = LensletGridModel.VSpacing + mod(LensletGridModel.VSpacing,2); % make sure it is even
else
    LensletGridModel.VSpacing = LensletGridModel.HSpacing;
end
LensletGridModel.HOffset = Camera.horizOffset * superResFactor;
LensletGridModel.VOffset = Camera.vertOffset * superResFactor;
LensletGridModel.Rot = Camera.gridRot;