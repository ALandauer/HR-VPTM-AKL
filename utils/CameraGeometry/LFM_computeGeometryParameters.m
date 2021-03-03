% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function [LensletCenters, Resolution, LensletGridModel, NewLensletGridModel] = ...
                                                    LFM_computeGeometryParameters(Camera, WhiteImage, depthRange, depthStep, superResFactor, DebugBuildGridModel, imgSize)
%% Process white image to find the real lenslet centers
if isempty(WhiteImage) % for simulation purposes; build LensetGridModel from specs and imgSize (when a white image does not exist)
    
    % Build lenslet grid model (MLA descriptor)
    LensletGridModel = struct;
    MLASize = ceil(imgSize/Camera.newSpacingPx); % effective size of the microlens array (in lenslets)
    LensletGridModel.UMax = MLASize(2); % no of lenslets
    LensletGridModel.VMax = MLASize(1);

    if(strcmp(Camera.gridType, 'hex'))
        LensletGridModel.VSpacing = round(sqrt(3)/2*Camera.spacingPx);
    else
        LensletGridModel.VSpacing = round(Camera.spacingPx);
    end
    LensletGridModel.HSpacing = round(Camera.spacingPx);
    LensletGridModel.FirstPosShiftRow = 1;
    LensletGridModel.Orientation = 'horz';
    LensletGridModel.HOffset = 0; %mod(ceil(imgSize(2)/2), Camera.newSpacingPx);
    LensletGridModel.VOffset = 0; %mod(ceil(imgSize(1)/2), Camera.newSpacingPx);
else
    LensletGridModel = LFM_processWhiteImage(WhiteImage, Camera.spacingPx, Camera.gridType, DebugBuildGridModel);
end

%% Transform to integer centers position 
% create the desired grid model when choosing a new spacing between
% lenslets (in pixels); Camera.newSpacingPx determines the desired sensor resolution
HOffset = 0; VOffset = 0; Rot = 0;
[NewLensletGridModel] = LFM_setGridModel(Camera.newSpacingPx, LensletGridModel.FirstPosShiftRow, LensletGridModel.UMax, LensletGridModel.VMax,...
                            HOffset, VOffset, Rot, LensletGridModel.Orientation, Camera.gridType);
                        
InputSpacing = [LensletGridModel.HSpacing, LensletGridModel.VSpacing];
NewSpacing = [NewLensletGridModel.HSpacing, NewLensletGridModel.VSpacing];
XformScale = NewSpacing ./ InputSpacing;  % Notice the resized image will not be square
NewOffset = round([LensletGridModel.HOffset, LensletGridModel.VOffset] .* XformScale);
NewLensletGridModel.HOffset = NewOffset(1);
NewLensletGridModel.VOffset = NewOffset(2);

%% Compute the texture side grid model 
if strcmp(superResFactor, 'default')
    superResFactor = Camera.newSpacingPx;
end
% When superResFactor == 1, we reconstruct at lenslet resolution
% When superResFactor == Nnum, we reconstruct at sensor resolution
[TextureGridModel] = LFM_setGridModel(superResFactor, LensletGridModel.FirstPosShiftRow, LensletGridModel.UMax, LensletGridModel.VMax,...
                            HOffset, VOffset, Rot, LensletGridModel.Orientation, Camera.gridType);

%% Compute resolution according to the new grid
Resolution = LFM_computeResolution(NewLensletGridModel, TextureGridModel, Camera, depthRange, depthStep);
Resolution.superResFactor = superResFactor;

disp(['Super resolution factor of: ', num2str(Resolution.TexNnum),' Pix size: [', num2str(Resolution.sensorRes(1)),', ',num2str(Resolution.sensorRes(2)),...
        '] Vox size: [', num2str(Resolution.texRes(1)), ', ',num2str(Resolution.texRes(2)),', ',num2str(Resolution.texRes(3)),']']);
    
%% Compute lenslets centers on the sensor and corresponding repetition patches centers in texture space (in voxels)

NewLensletGridModel.FirstPosShiftRow = LensletGridModel.FirstPosShiftRow;
TextureGridModel.FirstPosShiftRow = NewLensletGridModel.FirstPosShiftRow;

LensletCenters = LFM_computeLensCenters(NewLensletGridModel, TextureGridModel, Resolution.sensorRes, Camera.focus, Camera.gridType); 