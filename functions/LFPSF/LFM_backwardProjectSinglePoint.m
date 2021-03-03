% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function Backprojection = LFM_backwardProjectSinglePoint(H, Resolution, imgSize, texSize, currentPixel, lensletCenters, range, lensOrder)
% backwardProjectSinglePoint: This function computes the object response to a single senxsor pixel (which voxels from the object affect a single pixel in the sensor). 
% It return a stack of images equal to the number of depths.

if (ismatrix(H) && size(H,2) == 3) % for multi focus microscope
    multiFocus = 3;
    nDepths = size(H{1,1,1},3);
    ulensType = lensletCenters.px(:,:,3);
else
    multiFocus = 1;
    nDepths = size(H,3);
end

Backprojection = zeros(texSize(1), texSize(2), nDepths);
% Fetch lenslet centers in object space
lensVoxY = lensletCenters.vox(:,:,1);
lensVoxX = lensletCenters.vox(:,:,2);

% Fetch lenslet centers in image space
lensSenY = lensletCenters.px(:,:,1);
lensSenX = lensletCenters.px(:,:,2);

% Iterate different types of micro-lenses
for i = 1 : multiFocus
    if(multiFocus == 1)
        currentLensVoxY = lensVoxY;
        currentLensVoxX = lensVoxX;
        
        currentLensSenY = lensSenY;
        currentLensSenX = lensSenX;
    else
        currentLensVoxY = lensVoxY(ulensType == i);
        currentLensVoxX = lensVoxX(ulensType == i);
        
        currentLensSenY = lensSenY(ulensType == i);
        currentLensSenX = lensSenX(ulensType == i);
    end
    
    % Iterate all depths
    for cc = 1:nDepths
        sliceCurrentDepth = zeros(texSize);
        % Iterate voxels behind central lenslet
        for aa = 1:Resolution.TexNnum(1)
            aa_new = aa;
            flipX = 0;
            if (aa > Resolution.TexNnum_half(1) && strcmp(range, 'quarter'))
                aa_new = Resolution.TexNnum(1) - aa + 1;
                flipX = 1;
            end
            for bb = 1:Resolution.TexNnum(2)
                
                % Avoid overlaps for non-regular grids
                if  Resolution.texMask(aa,bb) == 0
                    continue;
                end
                
                bb_new = bb;
                flipY = 0;
                if (bb > Resolution.TexNnum_half(2) && strcmp(range, 'quarter'))
                    bb_new = Resolution.TexNnum(2) - bb + 1;
                    flipY = 1;
                end
                
                % Rotate forward pattern for convolution
                if(multiFocus == 1)
                    Ht = imrotate(H{aa_new,bb_new,cc}, 180);
                else
                    Ht = imrotate(H{1,lensOrder(i)}{aa_new,bb_new,cc}, 180);
                end
                
                % Convolve rotated forward pattern with a single poing, in image space
                tempSlice = sconv2singlePointFlip(imgSize, currentPixel, Ht, flipX, flipY, 'same');
                
                % Grab relevant voxels
                lensYInsideTex = round(currentLensVoxY - Resolution.TexNnum_half(1) + aa);
                lensXInsideTex = round(currentLensVoxX - Resolution.TexNnum_half(2) + bb);
                
                validLensTex = (lensXInsideTex <= texSize(2)) & (lensXInsideTex > 0) & ...
                    (lensYInsideTex <= texSize(1)) & (lensYInsideTex > 0);
                
                % Method 2 without resizing, 2x faster
                lensYInside = currentLensSenY + round((-Resolution.TexNnum_half(1) + aa) / Resolution.texScaleFactor(1));
                lensXInside = currentLensSenX + round((-Resolution.TexNnum_half(2) + bb) / Resolution.texScaleFactor(2));
                
                validLens = (lensXInside <= imgSize(2)) & (lensXInside > 0) & ...
                    (lensYInside <= imgSize(1)) & (lensYInside > 0);
                
                validLens = validLens & validLensTex;
                validLensTex = validLensTex & validLens;
                
                indicesSen = sub2ind(imgSize, lensYInside(validLens), lensXInside(validLens));
                indicesTex = sub2ind(texSize, lensYInsideTex(validLensTex), lensXInsideTex(validLensTex));
                
                if length(indicesSen) ~= length(indicesTex)
                    disp('Error in backProjectSinglePoint: mismatching indices!');
                end
                
                sliceCurrentDepth(indicesTex) = sliceCurrentDepth(indicesTex) + tempSlice(indicesSen);
            end
        end
        Backprojection(:,:,cc) = Backprojection(:,:,cc) + sliceCurrentDepth;
    end
end
