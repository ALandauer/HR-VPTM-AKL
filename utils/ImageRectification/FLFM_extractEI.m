% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu

function LF = FLFM_extractEI(LensletGridModel, LensletImage)
% Adapted after:
% LF Toolbox v0.4 released 12-Feb-2015
% Copyright (c) 2013-2015 Donald G. Dansereau

fprintf('\nSlicing lenslet image into LF elemental (sub-aperture) images.');

USize = LensletGridModel.UMax;
VSize = LensletGridModel.VMax;
MaxSpacing = max(LensletGridModel.HSpacing, LensletGridModel.VSpacing);  % Enforce square output in s,t
SSize = ceil(MaxSpacing/2)*2; %MaxSpacing - 1 - mod(MaxSpacing, 2); % force odd for centered middle pixel -- H,VSpacing are even, so +1 is odd
TSize = ceil(MaxSpacing/2)*2; %MaxSpacing - 1 - mod(MaxSpacing, 2);

SSize = SSize - 1 + mod(SSize, 2);
TSize = TSize - 1 + mod(TSize, 2);
LF = zeros(TSize, SSize, VSize, USize);

TVec = cast(floor((-(TSize-1)/2):((TSize-1)/2)), 'int16');
SVec = cast(floor((-(SSize-1)/2):((SSize-1)/2)), 'int16');
VVec = cast(0:VSize-1, 'int16');
UBlkSize = 32;
for( UStart = 0:UBlkSize:USize-1 ) 
    UStop = UStart + UBlkSize - 1;
    UStop = min(UStop, USize-1);  
    UVec = cast(UStart:UStop, 'int16');
    
    [tt,ss,vv,uu] = ndgrid( TVec, SVec, VVec, UVec );
    
    %---Build indices into 2D image---
    LFSliceIdxX = LensletGridModel.HOffset + uu.*LensletGridModel.HSpacing + ss;
    LFSliceIdxY = LensletGridModel.VOffset + vv.*LensletGridModel.VSpacing + tt;
    
    HexShiftStart = LensletGridModel.FirstPosShiftRow;
    
    %% TODO: address reg grid
    if (strcmp(LensletGridModel.gridType, 'hex'))
        LFSliceIdxX(:,:,HexShiftStart:2:end,:) = LFSliceIdxX(:,:,HexShiftStart:2:end,:) + LensletGridModel.HSpacing/2;
    end
    
    %---Lenslet mask in s,t and clip at image edges---
    R = sqrt(double(tt).^2 + double(ss).^2);
    ValidIdx = find(R < LensletGridModel.HSpacing/2 & ...
        LFSliceIdxX >= 1 & LFSliceIdxY >= 1 & LFSliceIdxX <= size(LensletImage,2) & LFSliceIdxY <= size(LensletImage,1) );
    
    %--clip -- the interp'd values get ignored via ValidIdx--
    LFSliceIdxX = max(1, min(size(LensletImage,2), LFSliceIdxX ));
    LFSliceIdxY = max(1, min(size(LensletImage,1), LFSliceIdxY ));

    LFSliceIdx = sub2ind(size(LensletImage), cast(LFSliceIdxY,'int32'), ...
        cast(LFSliceIdxX,'int32'), ones(size(LFSliceIdxX),'int32'));
    
    tt = tt - min(tt(:)) + 1;
    ss = ss - min(ss(:)) + 1;
    vv = vv - min(vv(:)) + 1;
    uu = uu - min(uu(:)) + 1 + UStart;
    LFOutSliceIdx = sub2ind(size(LF), cast(tt,'int32'), cast(ss,'int32'), ...
        cast(vv,'int32'),cast(uu,'int32'), ones(size(ss),'int32'));

        LF(LFOutSliceIdx(ValidIdx)) = ...
            LensletImage( LFSliceIdx(ValidIdx));
    fprintf('.');
end
end
