% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function [LensletImage, WhiteImage, configFile] = LFM_selectImages(dataset)

%% Choose the raw image to reconstruct and correponding white image (to identify the mlens centers)
%folder = '../../SampleData/LFM/';

folder = 'S:\Individual\Selda\LFM\LightFieldImaging_FinalPaperWork\Olaf\olaf-master\olaf-master\SampleData\LFM\';

subfolder = dataset;
% files = dir(fullfile([folder, subfolder,'/'], '*.tif'));

LensletImage = imread([folder, subfolder,  '\LFImage.tif']);
WhiteImage = imread([folder, subfolder,  '\WhiteImage.tif']);
configFile = [folder, subfolder,  '\LFMconfig.yaml'];

% figure; imagesc(LensletImage);
% rect = round(getrect)

% crop files to ROI before reconstructing for speed up
if strcmp(dataset, 'fishEye')
    rect = [1154, 270, 1200, 1227];
elseif strcmp(dataset, 'organoid')
    rect = [20, 25, 471, 461];
elseif strcmp(dataset, 'dx_exp')
    rect = [306,117,717,545];
elseif strcmp(dataset, 'fishMulti')
    rect = [614, 1132, 449, 496];   
elseif strcmp(dataset, 'synthetic')
    rect = [0, 0, 1280, 800];   
else % spheres
%     rect = [374 359 1261 1367];
    rect = [971, 363, 707, 1026];
end

LensletImage = imcrop(LensletImage, rect);
WhiteImage = imcrop(WhiteImage, rect);