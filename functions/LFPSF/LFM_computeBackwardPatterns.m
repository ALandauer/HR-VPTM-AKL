% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function Ht = LFM_computeBackwardPatterns(H, Resolution, range, lensOrder)
% ComputeBackwardPatterns: Computes which light points in the
% object affect every pixel behind a micro-lens.

%% retrieve sensor image and 3D scene containers sizes
if (ismatrix(H) && size(H,2) == 3) % multi focus MLA case 
    nDepths = size(H{1,1},3); % depths
    imgSize = size(H{1,1}{Resolution.TexNnum_half(1), Resolution.TexNnum_half(2)}); % forward projection size
else
    nDepths = size(H,3); % depths
    imgSize = size(H{1,1,1}); % forward projection size
end

% Compute volume size; it will be different to the image size in case of a
% superRFactor different than the number of pixels behind a lenslet.
texSize = ceil(imgSize.*Resolution.texScaleFactor);
texSize = texSize + (1-mod(texSize,2));

%% offset the lenslet centers to match the image/volume centers
offsetImg = ceil(imgSize./2); 
offsetVol = ceil(texSize./2);
lensletCenters.px(:,:,1) = Resolution.usedLensletCenters.px(:,:,1) + offsetImg(1);
lensletCenters.px(:,:,2) = Resolution.usedLensletCenters.px(:,:,2) + offsetImg(2);
if (ismatrix(H) && size(H,2) == 3)
    lensletCenters.px(:,:,3) = Resolution.usedLensletCenters.px(:,:,3);
end
lensletCenters.vox(:,:,1) = Resolution.usedLensletCenters.vox(:,:,1) + offsetVol(1);
lensletCenters.vox(:,:,2) = Resolution.usedLensletCenters.vox(:,:,2) + offsetVol(2);

%% for regular grids compute the backprojection patterns only for one quarter of coordinates (due to symmetry)
if strcmp(range, 'quarter')
    coordsRange  = Resolution.Nnum_half;
else
    coordsRange  = Resolution.Nnum;
end

%% compute backprojection patterns for all the pixels in coordsRange
Ht = cell(coordsRange(1), coordsRange(2), nDepths); % container for back projection patterns

% Iterate through every pixel behind the central micro-lens and compute
% which part of the texture (object) affects it.
for aa_sensor = 1:coordsRange(1)
    aa_tex = ceil(aa_sensor * Resolution.texScaleFactor(1)); % compute the corresponding coord of "aa_sensor" pixel in real world space
    parfor bb_sensor = 1:coordsRange(2)
        bb_tex = ceil(bb_sensor * Resolution.texScaleFactor(2)); % compute the corresponding coord of "bb_sensor" pixel in real world space
        
        % backproject the activated sensor pixel 
        currentPixel = [aa_sensor + offsetImg(1) - Resolution.Nnum_half(1), bb_sensor + offsetImg(2) - Resolution.Nnum_half(2)]; % position of the current active pixel in the image
        tempback = LFM_backwardProjectSinglePoint(H, Resolution, imgSize, texSize, currentPixel, lensletCenters, range, lensOrder);
        
        % bring the backprojection to center (the operator assumes the bproj patterns are centered when applying them)
        tempbackShift = cell(1, nDepths);
        for cc = 1:nDepths
            backShiftX = round(Resolution.TexNnum_half(1)-aa_tex);
            backShiftY = round(Resolution.TexNnum_half(2)-bb_tex); 
            shifted = imShift2(tempback(:,:,cc), backShiftX, backShiftY);
            tempbackShift{1,cc} = sparse(shifted);
        end
        
        % store the pattern
        Ht(aa_sensor,bb_sensor,:) = tempbackShift;
    end
    disp(['sensorX: ',num2str(aa_sensor)]);
end
