% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

% LFM_backwardProject: back projects a lenslet image into a volume.
function BackProjection = LFM_backwardProject(Ht, projection, lensCenters, Resolution, texSize, range)

Nnum = Resolution.Nnum;
Nnum_half = Resolution.Nnum_half;

%% retrieve sensor image and 3D scene containers sizes
nDepths = size(Ht,3);
imgSize = size(projection);

% offset centers to match the image and 3D backprojection
offsetImg = ceil(imgSize./2); 
offsetVol = ceil(texSize./2);
lensXpx = lensCenters.px(:,:,1) + offsetImg(1);
lensYpx = lensCenters.px(:,:,2) + offsetImg(2);

%% precompute positions in the vol/image where to apply the aa/bb backprojection patterns to/from
zeroSlice = zeros(texSize);
indicesImg = cell(Nnum(1),Nnum(2));
indicesTex = cell(Nnum(1),Nnum(2));
tempSlices = cell(Nnum(1),Nnum(2));

for aa_sen = 1:Nnum(1)
    for bb_sen = 1:Nnum(2)
        lensXpxCurrent = round(lensXpx - Nnum_half(1) + aa_sen); % why x becomes y here?
        lensYpxCurrent = round(lensYpx - Nnum_half(2) + bb_sen);

        % Find corresponding voxel for a given pixel. In case of different
        % voxel/pixel number more than one pixel will hit a voxel
        lensXvoxCurrent = lensXpxCurrent - offsetImg(1);
        lensXvoxCurrent = lensXvoxCurrent*Resolution.texScaleFactor(1);
        lensXvoxCurrent = ceil(lensXvoxCurrent + offsetVol(1));
        
        lensYvoxCurrent = lensYpxCurrent - offsetImg(2);
        lensYvoxCurrent = lensYvoxCurrent*Resolution.texScaleFactor(2);
        lensYvoxCurrent = ceil(lensYvoxCurrent + offsetVol(2));

        % check for out of image and texture
        validLens = (lensXpxCurrent <= imgSize(1)) & (lensXpxCurrent > 0) & ...
            (lensYpxCurrent <= imgSize(2)) & (lensYpxCurrent > 0) & ...
            (lensXvoxCurrent <= texSize(1)) & (lensXvoxCurrent > 0) & ...
            (lensYvoxCurrent <= texSize(2)) & (lensYvoxCurrent > 0);
        
        lensXpxCurrent = lensXpxCurrent(validLens);
        lensYpxCurrent = lensYpxCurrent(validLens);
        indicesImg{aa_sen,bb_sen} = sub2ind(imgSize, lensXpxCurrent, lensYpxCurrent);
              
        lensXvoxCurrent = lensXvoxCurrent(validLens);
        lensYvoxCurrent = lensYvoxCurrent(validLens);
        indicesTex{aa_sen,bb_sen} = sub2ind(texSize, lensXvoxCurrent, lensYvoxCurrent);
        
        % Sample image and store the corresponing pixels to project
        tempSlice = zeroSlice;
        tempSlice(indicesTex{aa_sen,bb_sen}) = projection(indicesImg{aa_sen,bb_sen});
        tempSlices{aa_sen,bb_sen} = tempSlice;
    end
end

%% Backproject
BackProjection = zeros(texSize(1),texSize(2), nDepths);

if sum(projection(:)) == 0
    return
end

% iterate all depths
for cc = 1:nDepths
    fprintf(['\nBP depth: ', num2str(cc), '/', num2str(nDepths)]);
    tempSliceBack = zeroSlice;
    
    for aa_sen = 1:Nnum(1)
        
        % precompute pattern index (for regular grids compute the backprojection patterns only for one quarter of coordinates (due to symmetry))
        aa_new = aa_sen;
        flipX = 0;
        if (strcmp(range, 'quarter') && aa_sen > Nnum_half(1))
            aa_new = Nnum(1) - aa_sen + 1;
            flipX = 1;
        end
        
        % fetch only necesary PSF's to reduce the amount copied to the
        % threads in the parfor
        HtcurrentAACC = Ht(aa_new,:,cc);
        tempSlicesToUse = tempSlices(aa_sen,:);
        
        parfor bb_sen = 1:Nnum(2)    
            % Backproject every pixel once. Avoid overlap (hex grid)
            if Resolution.sensMask(aa_sen, bb_sen) == 0
                continue;
            end
            
            % precompute pattern index (for regular grids compute the backprojection patterns only for one quarter of coordinates (due to symmetry))
            bb_new = bb_sen;
            flipY = 0;
            if (strcmp(range, 'quarter') && bb_sen > Nnum_half(2))
                bb_new = Nnum(2) - bb_sen + 1;
                flipY = 1;
            end
               
            % Instead of iterating through every lenslet, manage them as a
            % whole, Matlab is more efficient this way
            Hts = HtcurrentAACC{bb_new};
            if flipX ~= 0
                Hts = flipud(Hts);
            end
            if flipY ~= 0
                Hts = fliplr(Hts);
            end
            
            % Fetch sampled image at aa/bb pixels
            tempSlice = tempSlicesToUse{bb_sen};
            
            % convolve
            if sum(tempSlice(:)) > 0
                tempSliceBack = tempSliceBack + conv2(tempSlice, full(Hts), 'same');
            end
        end
    end   
        BackProjection(:,:,cc) = BackProjection(:,:,cc) + tempSliceBack;
end