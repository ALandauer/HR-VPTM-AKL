% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

% Based on:
% LF Toolbox v0.4 released 12-Feb-2015
% Copyright (c) 2013-2015 Donald G. Dansereau
function [CorrectedLensletImage, CorrectedWhiteImage] = LFM_applyTransformation(LensletImage, WhiteImage, FixAll, LensletCenters, debug)

XformScale = [FixAll.tdata.T(1,1), FixAll.tdata.T(2,2)];
NewSize = round(size(LensletImage(:,:,1)) .* XformScale(2:-1:1));
NewSize = floor(NewSize/2)*2 + 1;
CorrectedLensletImage = imtransform( LensletImage, FixAll, 'YData',[1 NewSize(1)], 'XData',[1 NewSize(2)]);
CorrectedWhiteImage = imtransform( WhiteImage, FixAll, 'YData',[1 NewSize(1)], 'XData',[1 NewSize(2)]);

% account for the center offset 
offset = ceil(size(CorrectedWhiteImage)/2) - LensletCenters.offset;
CorrectedWhiteImage = imtranslate(CorrectedWhiteImage, [offset(2), offset(1)], 'linear');
CorrectedLensletImage = imtranslate(CorrectedLensletImage, [offset(2), offset(1)], 'linear');

if (debug)
    lensX = LensletCenters.px(:,:,1) + size(CorrectedLensletImage,1)/2;
    lensY = LensletCenters.px(:,:,2) + size(CorrectedLensletImage,2)/2;

    figure; imagesc(CorrectedWhiteImage);
    hold on;
    scatter(lensY(:),lensX(:))
    hold off;
    title ('rectified centers');
    drawnow
end