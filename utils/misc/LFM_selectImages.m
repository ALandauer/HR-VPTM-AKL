% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page
%
% Upated for image sequence processing for LFM2PT, AKL 03/2021
%

function [LensletImageSeq, image_names, WhiteImage, configFile] = LFM_selectImages(folder,subfolder)

% Choose the raw images to reconstruct and correponding white image (to identify the mlens centers)
files = dir(fullfile(folder, subfolder, '*.tif'));

cnt = 1;
for ii = 1:length(files)
    if ~strcmpi(files(ii).name,'WhiteImage.tif')
        LensletImageSeq{cnt} = imread(fullfile(files(ii).folder,files(ii).name));
        image_names{cnt} = fullfile(files(ii).folder,files(ii).name);
        cnt = cnt + 1;
    end
end
WhiteImage = imread([folder, subfolder,  '\WhiteImage.tif']);
configFile = [folder, subfolder,  '\LFMconfig.yaml'];

%show figure for cropping
[~,rect] = imcrop(LensletImageSeq{1}+WhiteImage/4);
[LensletImageSeq{1}] = imcrop(LensletImageSeq{1},rect);

% crop files to ROI before reconstructing for speed up
if length(length(LensletImageSeq))>1
    for ii = 2:length(LensletImageSeq)
        LensletImageSeq{ii} = imcrop(LensletImageSeq{ii}, rect);
    end
end
WhiteImage = imcrop(WhiteImage, rect);

close(gcf)
