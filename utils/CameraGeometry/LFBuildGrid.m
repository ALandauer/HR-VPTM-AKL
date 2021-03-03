% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

% Adapted (to accomodate for regular grid arrays) after:
% LF Toolbox v0.4 released 12-Feb-2015
% Copyright (c) 2013-2015 Donald G. Dansereau

function [GridCoords] = LFBuildGrid(LensletGridModel, gridType)

RotCent = eye(3);
RotCent(1:2,3) = [LensletGridModel.HOffset, LensletGridModel.VOffset];

ToOffset = eye(3);
ToOffset(1:2,3) = [LensletGridModel.HOffset, LensletGridModel.VOffset];

R = ToOffset * RotCent * eul2rotm([LensletGridModel.Rot, 0, 0]) * RotCent^-1;

[vv,uu] = ndgrid((0:LensletGridModel.VMax-1).*LensletGridModel.VSpacing, (0:LensletGridModel.UMax-1).*LensletGridModel.HSpacing);

if(strcmp(gridType, 'hex'))
    uu(LensletGridModel.FirstPosShiftRow:2:end,:) = uu(LensletGridModel.FirstPosShiftRow:2:end,:) + 0.5.*LensletGridModel.HSpacing;
end

GridCoords = [uu(:), vv(:), ones(numel(vv),1)];
GridCoords = (R*GridCoords')';
GridCoords = reshape(GridCoords(:,1:2), [LensletGridModel.VMax,LensletGridModel.UMax,2]);
