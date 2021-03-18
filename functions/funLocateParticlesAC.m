function [x, beadParameter] = funLocateParticlesAC(I, beadParameter,vol_num)
% [x, beadParameter] = funLocateParticlesAC(I, beadParameter) locates particles in the image
%
% INPUTS
% -------------------------------------------------------------------------
%   I:              Input volumetric image
%   beadParameter:  Parameters to detect particle position in images (see 'funSetUpBeadParams.m')
%   vol_num:        Number of volume in the sequence
%
% OUTPUTS
% -------------------------------------------------------------------------
%   x:              Voxel-level estimate of particle center in MxNxO format
%

% Parameters
thres = beadParameter.thres;    %Threshold value
ratThresh = beadParameter.ratThresh; %axial to in-plane threshhold
smoothFac = beadParameter.smoothFac; %active contour smoothing factor
circThresh = beadParameter.circThresh;

% other params, set by user set later in this method
% minPixels = beadParameter.minSize;  %Minimum pixel count in blob for bead
% maxPixels = beadParameter.maxSize;  %Maximum pixel count in blob for bead

% %normalize
% I = I/max(I(:));% commented, should already be normalized coming from Trial-MPT

disp('%%%%%% Starting Binarization %%%%%%')

% Image thresholding
vol_binary = imbinarize(I,thres); %threshold-based to start

%dilate the result with a relatively large structuring element to
%privide a good initial guess to active contours
se_bin = strel('sphere',4);
% se8 = strel('sphere',8);
L = imdilate(vol_binary, se_bin); %if using watershed, replace "vol_binary" with "L"


%use active contour segmentation to refine binarization
BWseg = activecontour(I,L,100,'Chan-Vese','SmoothFactor',smoothFac);
se_seg = strel('sphere',1);
BWseg = imdilate(BWseg, se_seg); %dilate slightly to improve connectivity

% Find bead blobs
CC = bwconncomp(BWseg);
numPixels = cellfun(@numel,CC.PixelIdxList);

if vol_num == 1
    %check the histogram of sizes and get user defined size limits
    try
    nbins = sshist(numPixels);
    catch
        nbins = 1;
    end
    nbins = max([15,nbins]);
    
    h = figure;
    histogram(numPixels,nbins)
    minPixels = input('Enter min bead size: ');  %Minimum pixel count in blob for bead
    maxPixels = input('Enter max bead size: ');  %Maximum pixel count in blob for bead
    try
        close(h)
    catch
    end
% %     minPixels = 70;
% %     maxPixels = 1250;
% %     beadParameter.minSize = minPixels;  %Minimum pixel count in blob for bead
% %     beadParameter.maxSize = maxPixels;  %Maximum pixel count in blob for bead
    beadParameter.minSize = minPixels;  %Minimum pixel count in blob for bead
    beadParameter.maxSize = maxPixels;  %Maximum pixel count in blob for bead
    
else
    
    minPixels = beadParameter.minSize;  %Minimum pixel count in blob for bead
    maxPixels = beadParameter.maxSize;  %Maximum pixel count in blob for bead
    
end
% Remove beads from the list that are bigger and smaller than the limits
beadBlob = numPixels>minPixels & numPixels<maxPixels;

% Find centroid of all the eligible blobs
% get region info
S_ = regionprops3(CC, 'Centroid', 'EigenValues', 'PrincipalAxisLength', 'VoxelList');

%calculate mean essentricity
try
    AxisLenRatXZ = S_.PrincipalAxisLength(:,2)./S_.PrincipalAxisLength(:,1);
    AxisLenRatYZ = S_.PrincipalAxisLength(:,3)./S_.PrincipalAxisLength(:,1);
    circularity = S_.PrincipalAxisLength(:,2)./S_.PrincipalAxisLength(:,3) - 1;
    mean_rat_Z = mean([AxisLenRatXZ,AxisLenRatYZ],2);
catch
    mean_rat_Z = 0.7;
end

%get minimal region region props
s = regionprops(CC, 'Centroid');

%get centroid points of the eligible blobs
blobPts = round(double(struct2dataset(s)));
blobPts = blobPts(beadBlob' & mean_rat_Z>ratThresh & abs(circularity)<circThresh,:);
temp = blobPts;

% Convert to m,n,o coordinates
blobPts(:,1) = temp(:,2);
blobPts(:,2) = temp(:,1);
x = blobPts;

if size(blobPts,1) > 1
    disp(length(blobPts))
end

disp('%%%%%% Binarization complete! %%%%%%')

end

