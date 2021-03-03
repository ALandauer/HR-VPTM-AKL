% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function LensletGridModel = LFM_setGridModel(SpacingPx, FirstPosShiftRow, UMax, VMax, HOffset, VOffset, Rot, Orientation, gridType)

% todo: defaults
if(strcmp(gridType, 'hex'))
    Spacing = [SpacingPx*cosd(30), SpacingPx];
    Spacing = ceil(Spacing);
    Spacing = ceil(Spacing/2)*2;
end

if(strcmp(gridType, 'reg'))
    Spacing = [SpacingPx, SpacingPx];
end

LensletGridModel.HSpacing = Spacing(2);
LensletGridModel.VSpacing = Spacing(1);

LensletGridModel.HOffset = HOffset;
LensletGridModel.VOffset = VOffset;
LensletGridModel.Rot = Rot;
LensletGridModel.UMax = UMax;
LensletGridModel.VMax = VMax;
LensletGridModel.Orientation = Orientation;
LensletGridModel.FirstPosShiftRow = FirstPosShiftRow;