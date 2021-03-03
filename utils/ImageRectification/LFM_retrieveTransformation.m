% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

% Based on:
% LF Toolbox v0.4 released 12-Feb-2015
% Copyright (c) 2013-2015 Donald G. Dansereau

function FixAll = LFM_retrieveTransformation(LensletGridModel, NewLensletGridModel)

%% Notice vertical/horizontal inverted ->the toolbox uses a different conversion
% scale transform
InputSpacing = [LensletGridModel.HSpacing, LensletGridModel.VSpacing];
NewSpacing = [NewLensletGridModel.HSpacing, NewLensletGridModel.VSpacing];
XformScale = NewSpacing ./ InputSpacing;  % Notice the resized image will not be square

RScale = eye(3);
RScale(1,1) = XformScale(1);
RScale(2,2) = XformScale(2);

NewOffset = [LensletGridModel.HOffset, LensletGridModel.VOffset] .* XformScale;
RoundedOffset = round(NewOffset);
XformTrans =  RoundedOffset-NewOffset;

RTrans = eye(3);
RTrans(end,1:2) = XformTrans;

% rotation transform
RRot = eul2rotm([LensletGridModel.Rot, 0, 0]);

% final transformation
FixAll = maketform('affine', RRot*RScale*RTrans);