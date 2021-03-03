% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

% LFM_forwardProject: Forward projects a volume to a lenslet image simulating the behavior of the microscope.
function Projection = LFM_forwardProject( H, realSpace, lensCenters, Resolution, imgSize, range)

TexNnum = Resolution.TexNnum;
TexNnum_half = Resolution.TexNnum_half;

%% Retrieve sensor image and 3D scene containers sizes
nDepths = size(H,3);
texSize = [size(realSpace,1), size(realSpace,2)];

% offset centers to match the image and 3D backprojection
offsetImg = ceil(imgSize./2);
offsetVol = ceil(texSize./2);
lensYvox = lensCenters.vox(:,:,1) + offsetVol(1);
lensXvox = lensCenters.vox(:,:,2) + offsetVol(2);

%% Precompute positions in the vol/image where to apply the different PSFs patterns
indicesTex = cell(TexNnum(1),TexNnum(2));
indicesImg = cell(TexNnum(1),TexNnum(2));

for aa_tex = 1:TexNnum(1)
    for bb_tex = 1:TexNnum(2)
        % Find voxels to sample relative to the lenslet centers in the
        % volume
        lensXvoxCurrent = round(lensYvox - TexNnum_half(1) + aa_tex);
        lensYvoxCurrent = round(lensXvox - TexNnum_half(2) + bb_tex);
        
        % convert lenses to img space
        lensXpxCurrent = lensXvoxCurrent - offsetVol(1);
        lensXpxCurrent = lensXpxCurrent/Resolution.texScaleFactor(1);
        lensXpxCurrent = ceil(lensXpxCurrent + offsetImg(1));
        
        lensYpxCurrent = lensYvoxCurrent - offsetVol(2);
        lensYpxCurrent = lensYpxCurrent/Resolution.texScaleFactor(2);
        lensYpxCurrent = ceil(lensYpxCurrent + offsetImg(2));
        
        % check for out of image and texture
        validLens = (lensXvoxCurrent <= texSize(1)) & (lensXvoxCurrent > 0) & ...
            (lensYvoxCurrent <= texSize(2)) & (lensYvoxCurrent > 0) & ...
            (lensXpxCurrent <= imgSize(1)) & (lensXpxCurrent > 0) & ...
            (lensYpxCurrent <= imgSize(2)) & (lensYpxCurrent > 0);
        
        lensXvoxCurrent = lensXvoxCurrent(validLens);
        lensYvoxCurrent = lensYvoxCurrent(validLens);
        indicesTex{aa_tex,bb_tex} = sub2ind(texSize,lensXvoxCurrent, lensYvoxCurrent);
        
        lensXpxCurrent = lensXpxCurrent(validLens);
        lensYpxCurrent = lensYpxCurrent(validLens);
        indicesImg{aa_tex,bb_tex} = sub2ind(imgSize,lensXpxCurrent, lensYpxCurrent);
    end
end

% Forwardproject
zeroSpace = zeros(imgSize);
Projection = zeroSpace;

for cc = 1:nDepths
%     tic
    fprintf(['\nFP depth: ', num2str(cc), '/', num2str(nDepths)]);
    realspaceCurrentDepth = realSpace(:,:,cc);
    
    if sum(realspaceCurrentDepth(:)) == 0
        continue;
    end
    for aa_tex = 1:TexNnum(1)
        
        % Precompute pattern index (for regular grids we computed the fwdprojection patterns only for one quarter of coordinates (due to symmetry))
        aa_new = aa_tex;
        flipX = 0;
        if (strcmp(range, 'quarter') && aa_tex > TexNnum_half(1))
            aa_new = TexNnum(1) - aa_tex + 1;
            flipX = 1;
        end
        
        % Slice PSF's needed for parallel computing
        HcurrentAACC = H(aa_new,:,cc);
        
            parfor bb_tex = 1:TexNnum(2)
                
                % Forward project from every point once. Avoid overlap (hex grid)
                if Resolution.texMask(aa_tex, bb_tex) == 0
                    continue;
                end
                
                % Precompute pattern index (for regular grids compute the backprojection patterns only for one quarter of coordinates (due to symmetry))
                bb_new = bb_tex;
                flipY = 0;
                if (strcmp(range, 'quarter') && bb_tex > TexNnum_half(2))
                    bb_new = TexNnum(2) - bb_tex + 1;
                    flipY = 1;
                end
                
                % Fetch corresponding PSF for given coordinate behind the
                % lenslet
                Hs = HcurrentAACC{bb_new};
                if flipX ~= 0
                    Hs = flipud(Hs);
                end
                if flipY ~= 0
                    Hs = fliplr(Hs);
                end
                
                tempspace = zeros(imgSize);
                tempspace(indicesImg{aa_tex,bb_tex}) = realspaceCurrentDepth(indicesTex{aa_tex,bb_tex});
                
                % Accumulate result
                if sum(tempspace(:)) > 0
                    projectedPattern = conv2(tempspace, full(Hs), 'same');
                    Projection = Projection + projectedPattern;
                end
            end
    end
end